#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
    //initialized
    is_initialized_ = false;
    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3 - n_x_;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  // initial state vector
  x_ = VectorXd(5);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5;
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;
  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
    weights_=VectorXd(2*n_aug_+1);
    previous_timestamp_=0;
    NIS_laser_=0.0;
    NIS_radar_=0.0;
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
    if(!is_initialized_){
        x_<<1, 1, 1, 1, 1;

        P_<<1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1;
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
            double ro=meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double ro_dot = meas_package.raw_measurements_(2);

            double vx=ro_dot*cos(phi);
            double vy=ro_dot*sin(phi);

            x_(0)=ro*cos(phi);
            x_(1)=ro*sin(phi);
            x_(2)=4;//sqrt(pow(vx,2) + pow(vy,2));
            x_(3)=ro_dot*cos(phi);//0;
            x_(4)=ro_dot*sin(phi);
            
//            P_ << std_radr_*std_radr_, 0, 0, 0, 0,
//            0, std_radr_*std_radr_, 0, 0, 0,
//            0, 0, 1, 0, 0,
//            0, 0, 0, std_radphi_, 0,
//            0, 0, 0, 0, std_radphi_;
            
        }
        else if(meas_package.sensor_type_ == MeasurementPackage::LASER){

            x_(0)=meas_package.raw_measurements_(0);
            x_(1)=meas_package.raw_measurements_(1);
            x_(2)=4;
            x_(3)=0.5;
            x_(4)=0;
            
//            P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
//            0, std_laspy_*std_laspy_, 0, 0, 0,
//            0, 0, 1, 0, 0,
//            0, 0, 0, 1, 0,
//            0, 0, 0, 0, 1;
        }

        previous_timestamp_ = meas_package.timestamp_;
        is_initialized_ =true;
        return;
    }
    //start the pipeline


    // Prediction
    double df = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(df);

    //update
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
//        cout<<"UpdateRadar-------"<<endl;
        UpdateRadar(meas_package);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
  //      cout<<"UpdateLidar-------"<<endl;
        UpdateLidar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //first generate sigma points
    //start first sigma points
    lambda_=3-n_x_;
    MatrixXd Xsig = MatrixXd(n_x_, 2*n_x_+1);
    MatrixXd A = P_.llt().matrixL();

    VectorXd x_aug=VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    // x_aug 7
    x_aug.head(5)=x_;
    x_aug(5)=0;
    x_aug(6)=0;
    // P_aug 7,7
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_)=P_;
    //std_a_ m/s^2
    P_aug(5,5)=std_a_*std_a_;
    //std_yawdd_  in rad/s^2
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    
    Xsig.col(0)=x_;
    for(int i=0;i<n_x_;i++){
        Xsig.col(i+1)=x_ + sqrt(lambda_ + n_x_)*A.col(i);
        Xsig.col(i+1+n_x_)=x_ - sqrt(lambda_ + n_x_)* A.col(i);
    }
    //second augmentation


    //create square root matrix
    MatrixXd A_aug = P_aug.llt().matrixL();
    //7,15
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
    lambda_=3-n_aug_;
    Xsig_aug.col(0) =x_aug;
    for (int i=0; i<n_aug_; i++) {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)* A_aug.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A_aug.col(i);
    }
    //sigma point prediction
    for(int i=0;i<2*n_aug_+1;i++){
        double p_x=Xsig_aug(0,i);
        double p_y=Xsig_aug(1,i);
        double v=Xsig_aug(2,i);
        double yaw=Xsig_aug(3,i);
        double yawd=Xsig_aug(4,i);
        double nu_a=Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        double px_p, py_p;
        if(fabs(yawd)>0.001){
            px_p = p_x+ v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw));
            py_p = p_y+ v/yawd*(-cos(yaw+yawd*delta_t)+cos(yaw));
        }else{
            px_p=p_x + v*cos(yaw)*delta_t;
            py_p=p_y+v*sin(yaw)*delta_t;
        }
        double v_p=v;
        double yaw_p = yaw+yawd*delta_t;//?
        double yawd_p = yawd;
        //add noise
        px_p = px_p + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
        py_p=py_p+0.5*delta_t*delta_t*sin(yaw)*nu_a;

        v_p = v_p + nu_a*delta_t;
        yaw_p = yaw_p + 0.5*delta_t*delta_t*nu_yawdd;
        yawd_p = yawd_p+ delta_t*nu_yawdd;

        //insert to sigma ,5,15
        Xsig_pred_(0,i)=px_p;
        Xsig_pred_(1,i)=py_p;
        Xsig_pred_(2,i)=v_p;
        Xsig_pred_(3,i)=yaw_p;
        Xsig_pred_(4,i)=yawd_p;
    }
    //set weights
    double weight_0=lambda_/(lambda_ + n_aug_);
    weights_(0)=weight_0;
    for(int i=1;i<2*n_aug_+1;i++){
        double weight=0.5/(n_aug_ + lambda_);
        weights_(i)= weight;
    }
    P_.fill(0);

    //Predict Mean and Covariance
    x_.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
//        mean
        x_ = x_ + weights_(i)*Xsig_pred_.col(i);
    }

    for(int i=0;i<2*n_aug_+1;i++){
        //        covariance

        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
        P_ =P_ + weights_(i)*x_diff*x_diff.transpose();
    }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    //2 measurement dimensions, lidar can px,py
    //lidar has 2 measurements
    int n_z=2;
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);//2,15
    VectorXd z_pred= VectorXd(n_z);//2
    MatrixXd S = MatrixXd(n_z, n_z);//2,2
    S.fill(0.0);
    z_pred.fill(0.0);
    Zsig.fill(0.0);

    MatrixXd R_laser_ = MatrixXd(n_z,n_z);
    R_laser_.fill(0.0);
    //create matrix for cross correlation Tc;
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    for(int i=0;i<2*n_aug_+1;i++){
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);

        Zsig(0, i)=p_x;
        Zsig(1, i) = p_y;
        //fuck
        z_pred += weights_(i) *Zsig.col(i);
    }
//    R_laser_(0,0)=std_laspx_*std_laspx_;
//    R_laser_(1,1) = std_laspy_*std_laspy_;
    R_laser_<<std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;

    for(int i=0;i<2*n_aug_+1;i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        while(z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
        while(z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;
        
        S = S+ weights_(i)*z_diff*z_diff.transpose();
        VectorXd x_diff = Xsig_pred_.col(i)-x_;
        while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
        while(x_diff(3)< -M_PI) x_diff(3) += 2.*M_PI;
        Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
    }
    S = S + R_laser_;
    //kalman gain K;
    MatrixXd K=Tc* S.inverse();
    VectorXd z_diff_ = meas_package.raw_measurements_ - z_pred;

    x_ =x_ + K * z_diff_;
    P_ = P_ - K * S * K.transpose();
    //Calculate NIS
    NIS_laser_ = z_diff_.transpose() * S.inverse() * z_diff_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    //3 measurement dimensions, radar can measure r, phi, and r_dot
    //radar has 3 measurements
    int n_z=3;
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
    VectorXd z_pred= VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    z_pred.fill(0.0);
    //UKF update assignment
    //create matrix for cross correlation Tc;
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    //init
    Tc.fill(0.0);

    //transform sigma points into measurement space
    for (int i = 0; i<2*n_aug_+1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);

        double vx=cos(yaw)*v;
        double vy=sin(yaw)*v;
        //rpo
        Zsig(0,i)=sqrt(px*px+py*py);
        Zsig(1,i)=atan2(py,px);
        Zsig(2,i)=(px*vx+py*vy)/sqrt(px*px+py*py);
        //
        z_pred = z_pred + weights_(i)*Zsig.col(i);
    }
    //1
    MatrixXd R = MatrixXd(n_z, n_z);//3,3

    R << std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;

    for(int i=0;i<2*n_aug_+1;i++){
        VectorXd z_diff = Zsig.col(i)-z_pred;
        while(z_diff(1) > M_PI) z_diff(1)-=2.*M_PI;
        while(z_diff(1) < -M_PI) z_diff(1)+=2.*M_PI;

        S = S+ weights_(i)*z_diff*z_diff.transpose();
        // to calculate the Tc(5,3)
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
        while(x_diff(3)<-M_PI) x_diff(3) +=2.*M_PI;

//        Tc += weights_(i)*x_diff*z_diff.transpose();
    }
    S += R;
    /* after calculate the S, update x & p
     */
    //calculate kalman gain K;

    MatrixXd K = Tc *S.inverse();
    VectorXd z = VectorXd(n_z);
    
    VectorXd z_diff_ = meas_package.raw_measurements_ - z_pred;
    while(z_diff_(1) > M_PI) z_diff_(1)-=2.*M_PI;
    while(z_diff_(1) < -M_PI) z_diff_(1)+=2.*M_PI;

    x_ = x_ + K * z_diff_;
    P_ = P_ - K * S * K.transpose();
    //Calculate NIS
    NIS_radar_ = z_diff_.transpose() * S.inverse() * z_diff_;
}
