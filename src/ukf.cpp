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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
    x_.fill(0.0);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);

    VectorXd x_aug(n_aug);
    MatrixXd P_aug = MatrixXd(n_aug, n_aug);
    
    x_aug.fill(0.0);
    P_aug.fill(0.0);
    
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
    P_<<1,0,0,0,0,
    0,1,0,0,0,
    0,0,1,0,0,
    0,0,0,1,0,
    0,0,0,0,1;
    
    int n_x=5;
    int n_aug=7;
    int lambda = 3- n_aug;
    VectorXd weights(2*n_aug+1);
    weights(0)=lambda/(lambda + n_aug);
    for(int i=1;i<n_aug;i++){
        weights(i)=0.5*(lambda+n_aug);
    }
    previous_timestamp_=0;
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2*n_aug+1);
    MatrixXd Xsig_pred = MatrixXd(n_x, 2*n_aug+1);
  
    //create vector for predicted state
    //create covariance matrix for prediction
    
    
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
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
            double ro=meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            double ro_dot = meas_package.raw_measurements_(2);
            
            double vx=ro_dot*cos(phi);
            double vy=ro_dot*sin(phi);
            x_(0)=ro*cos(phi);
            x_(1)=ro*sin(phi);
            x_(2)=sqrt(pow(vx,2) + pow(vy,2));
        }
        else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
            x_(0)=meas_package.raw_measurements_(0);
            x_(1)=meas_package.raw_measurements_(1);
            previous_timestamp_ = meas_package.timestamp_;
        }
        is_initialized_ =True;
        return;
    }
    //start the pipeline
    

    // Prediction
    double df = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
    
    Prediction(df);
    
    previous_timestamp_ = meas_package.timestamp_;
    
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
        UpdateRadar(meas_package);
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
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
    MatrixXd Xsig = MatrixXd(n_x, 2*n_x+1);
    MatrixXd A = P_.llt().matrixL();
    
    Xsig.col(0)=x_;
    for(int i=0;i<n_x;i++){
        Xsig.col(i+1)=x+sqrt(lambda + n_x)*A.col(i);
        Xsig.col(i+1+n_x)=x-sqrt(lambda + n_x)* A.col(i);
    }
    //second augmentation
    // x_aug 7
    x_aug.head(5)=x_;
    // P_aug 7,7
    P_aug.topLeftCorner(n_x, n_x)=P_;
    //std_a_ m/s^2
    P_aug(5,5)=std_a_*std_a_;
    //std_yawdd_  in rad/s^2
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    //create square root matrix
    MatrixXd Srm = P_aug.llt().matrixL();
    //7,15
    Xsig_aug.col(0) =x_aug;
    for (int i=0; i<n_aug; i++) {
        Xsig_aug.col(i+1) = x_aug+sqrt(lambda + n_aug)* Srm.col(i);
        Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda + n_aug)*Srm.col(i);
    }
    //sigma point prediction
    for(int i=0;i<2*n_aug+1;i++){
        double px=Xsig_aug(0,i);
        double py=Xsig_aug(1,i);
        double v=Xsig_aug(2,i);
        double yaw_angel=Xsig_aug(3,i);
        double yaw_rate=Xsig_aug(4,i);
        double nu_a=Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        double px_p, py_p;
        if(fabs(yaw_rate)>0.001){
            px_p = px+ v/yaw_rate*(sin(yaw_angel+yaw_rate*delta_t)-sin(yaw_angel));
            py_p = py+ v/yaw_rate*(-cos(yaw_angel+yaw_rate*delta_t)+cos(yaw_angel));
        }else{
            px_p=px + v*cos(yaw_angel)*delta_t;
            py_p=py+v*sin(yaw_angel)*delta_t;
        }
        double v_p=v;
        double yaw_angel_p = yaw_angel+yaw_angel*delta_t;
        double yaw_rate_p = yaw_rate;
        //add noise
        px_p = px_p + 0.5*delta_t*delta_t*cos(yaw_angel)*nu_a;
        py_p=py_p+0.5*delta_t*delta_t*sin(yaw_angel)*nu_a;
        
        v_p = v_p + nu_a*delta_t;
        yaw_angel_p = yaw_angel_p + 0.5*delta_t*delta_t*nu_yawdd;
        yaw_rate_p = yaw_rate_p+ delta_t*nu_yawdd;
        
        //insert to sigma ,5,15
        Xsig_pred(0,i)=px_p;
        Xsig_pred(1,i)=py_p;
        Xsig_pred(2,i)=v_p;
        Xsig_pred(3,i)=yaw_angel_p;
        Xsig_pred(4,i)=yaw_rate_p;
    }
    
    
    //Predict Mean and Covariance
    for(int i=0;i<2*n_aug+1;i++){
//        mean
        x_ = x_ + weights(i)*Xsig_pred.col(i);
    }
    for(int i=0;i<2*n_aug+1;i++){
        //        covariance
        VectorXd x_diff = Xsig_pred.col(i)-x_;
        while (x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
        P_ =P_ + weights(i)*x_diff*x_diff.transpose();
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
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug+1);//2,15
    VectorXd z_pred= VectorXd(n_z);//2
    MatrixXd S = MatrixXd(n_z, n_z);//2,2
    S.fill(0);
    z_pred.fill(0);
    Zsig.fill(0);
    
    MatrixXd R_laser_ = MatrixXd(n_z,n_z);
    R_laser_.fill(0.0);
    //create matrix for cross correlation Tc;
    MatrixXd Tc = MatrixXd(n_x, n_z);
    for(int i=0;i<2*n_aug+1;i++){
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        
        Zsig(0, i)=p_x;
        Zsig(1, i) = p_y;
    }
    R_laser_(0,0)=std_laspx_*std_laspx_;
    R_laser_(1,1) = std_laspy_*std_laspy_;
    
    for(int i=0;i<2*n_aug+1;i++){
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S+ weights(i)*z_diff*z_diff.transpose();
        
        VectorXd x_diff = Xsig_pred.col(i)-x_;
        Tc = Tc + weights(i)*x_diff*z_diff.transpose();
    }
    S += S;
    //kalman gain K;
    MatrixXd K=Tc* S.inverse();
    VectorXd z_diff_ = z - z_pred;
    
    x_ += K*z_diff_;
    P_ = P_ - K*S*K.transpose();
    
    //Calculate NIS
    double NIP = z_diff_.transpose() * S.inverse() * z_diff_;
    
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
    MatrixXd Zsig = MatrixXd(n_z, 2*n_aug+1);
    VectorXd z_pred= VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0);
    z_pred.fill(0);
    Zsig.fill(0);
    //UKF update assignment
    //create matrix for cross correlation Tc;
    MatrixXd Tc = MatrixXd(n_x, n_z);
    //init
    Tc.fill(0.0);
    
    //transform sigma points into measurement space
    for (int i = 0; i<2*n_aug+1; i++) {  //2n+1 simga points
        
        // extract values for better readibility
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v  = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);
        
        double vx=cos(yaw);
        double vy=sin(yaw);
        //rpo
        Zsig(0,i)=sqrt(px*px+py*py);
        Zsig(1,i)=atan2(py,px);
        Zsig(2,i)=(px*vx+py*vy)/sqrt(px*px+py*py);
    }
    MatrixXd R = MatrixXd(n_z, n_z);//3,3
    R.fill(0.0);
    R(0,0)=std_radr_*std_radr_;
    R(1,1)=std_radphi_*std_radphi_;
    R(2,2)=std_radrd_*std_radrd_;
    for(int i=0;i<2*n_aug+1;i++){
        z_pred += weights.col(i)*Zsig.col(i);
    }
    for(int i=0;i<2*n_aug+1;i++){
        VectorXd z_diff = Zsig.col(i)-z_pred;
        while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
        while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S+ weights(i)*z_diff*z_diff.transpose();
        // to calculate the Tc(5,3)
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
        while(x_diff(3)<-M_PI) x_diff(3) +=2.*M_PI;
        
        Tc = Tc + weights*x_diff*z_diff.transpose();
    }
    S += R;
    /* after calculate the S, update x & p
     */
    //calculate kalman gain K;
    MatrixXd K=Tc *S.inverse();
    VectorXd z_diff_=z-z_pred;
    while(z_diff_(1) > M_PI) z_diff_(1)-=2.*M_PI;
    while(z_diff_(1) < -M_PI) z_diff_(1)+=2.*M_PI;
    
    x=x+K*z_diff_;
    P_ = P_-K*S*K.transpose();
    
    //Calculate NIS
    double NIP = z_diff_.transpose() * S.inverse() * z_diff_;
}
