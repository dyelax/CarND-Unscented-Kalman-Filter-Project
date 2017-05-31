#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  n_x_ = 5;
  lambda_ = 3 - n_x_;

  n_aug_ = 7;
  n_sig_ = 2 * n_aug_ + 1;
  
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 16;

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
  
  // Init covar mat as identity (Might need to change some diag values?)
  P_ << 1, 0, 0,    0,    0,
        0, 1, 0,    0,    0,
        0, 0, 1000, 0,    0,
        0, 0, 0,    1,    0,
        0, 0, 0,    0,    1;
  
  // Predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  
  // Init weights
  weights_ = VectorXd(n_sig_);
  for (int i = 0; i < n_sig_; i++) {
    if (i == 0) {
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    } else {
      weights_(i) = 1 / (2 * (lambda_ + n_aug_));
    }
  }
}

UKF::~UKF() {}

VectorXd UKF::PolarToCartesian(const VectorXd &polar_measurements) {
  float rho = polar_measurements(0);
  float phi = polar_measurements(1);
  float rho_dot = polar_measurements(2);
  
  float px = rho * cos(phi);
  float py = rho * sin(phi);
//  float vx = rho_dot * cos(phi);
//  float vy = rho_dot * sin(phi);
  float v = rho_dot;
  
  VectorXd cartesian_measurements(3);
  cartesian_measurements << px, py, v;
  
  return cartesian_measurements;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  long long timestamp = meas_package.timestamp_;
  
  if (!is_initialized_) {
    previous_timestamp_ = timestamp;
    
    // first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      VectorXd cartesian_measurements =
      PolarToCartesian(meas_package.raw_measurements_);
      
      float px = cartesian_measurements(0);
      float py = cartesian_measurements(1);
      float v  = cartesian_measurements(2);
      
      x_ << px, py, v, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      x_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            0,
            0,
            0;
    }
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  // Get the time delta since the last measurement.
  double dt = (timestamp - previous_timestamp_) / 1000000.;  // Convert us to s.
  previous_timestamp_ = timestamp;
  
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
}

/** 
 * Computes the sigma points for an update step.
 */
MatrixXd UKF::GetSigPoints() {
  // Find sigma points.
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0);
  x_aug.head(5) = x_;
  // Mean value of noises (last 2 elts of augmented state) are 0, so leave them be.
  
  MatrixXd Q = MatrixXd(2, 2);
  Q << pow(std_a_, 2), 0,
       0,              pow(std_yawdd_, 2);
  
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q;
  
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  
  // Create square root matrix.
  MatrixXd A = P_aug.llt().matrixL();
  
  // Create augmented sigma points.
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1)         = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  return Xsig_aug;
}

/**
 * Applies the motion model f(x_k, nu_k) and returns x_(k+1).
 * @param x_aug an augmented state vector (usually a sigma point).
 * @param delta_t the time since the last prediction, in seconds
 */
VectorXd UKF::MotionModel(VectorXd &x_aug, double delta_t) {
  VectorXd covarTransVec = VectorXd(n_x_);
  VectorXd stateTransVec = VectorXd(n_x_);
  
  float px      = x_aug(0);
  float py      = x_aug(1);
  float v       = x_aug(2);
  float psi     = x_aug(3);
  float psid    = x_aug(4);
  float n_a     = x_aug(5);
  float n_psidd = x_aug(6);
  
  covarTransVec << 0.5 * pow(delta_t, 2) * cos(psi) * n_a,
                   0.5 * pow(delta_t, 2) * sin(psi) * n_a,
                   delta_t * n_a,
                   0.5 * pow(delta_t, 2) * n_psidd,
                   delta_t * n_psidd;
  
  float eps = 0.00001;
  if (fabs(psid) < eps) {
    stateTransVec << v * cos(psi) * delta_t,
                     v * sin(psi) * delta_t,
                     0,
                     psid * delta_t,
                     0;
  } else {
    stateTransVec << (v / psid) * (sin(psi + psid * delta_t) - sin(psi)),
                     (v / psid) * (-cos(psi + psid * delta_t) + cos(psi)),
                     0,
                     psid * delta_t,
                     0;
  }
  
  return x_aug.head(5) + stateTransVec + covarTransVec;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug = GetSigPoints();
  
  // Compute predictions on sigma points
  for (int i = 0; i < n_sig_; i++) {
    VectorXd sigpt = Xsig_aug.col(i);
    Xsig_pred_.col(i) = MotionModel(sigpt, delta_t);
  }

  // Predict mean and covar.
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) -= (2 * M_PI) * floor((x_diff(3) + M_PI) / (2 * M_PI));
    
    P_ += weights_(i) * x_diff * x_diff.transpose();
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
  
  // Measurement prediction for lidar

  // Update the state by using standard Kalman Filter equations.
  // Measurement model for lidar is linear wrt the state vector, so this is
  // equivalent to using sigma points, but cheaper.
  MatrixXd z = meas_package.raw_measurements_;
  int n_z = 2;
  
  MatrixXd H = MatrixXd(n_z, n_x_);
  MatrixXd R = MatrixXd(n_z, n_z);

  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  R << pow(std_laspx_, 2), 0,
       0,                  std_laspy_;
  
  
  // Intermediate calculations.
  MatrixXd Ht = H.transpose();
  VectorXd y = z - H * x_;
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();
  long size = x_.size();
  MatrixXd I = MatrixXd::Identity(size, size);
  
  // Update state and covariance mats.
  x_ = x_ + K * y;
  P_ = (I - K * H) * P_;
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
  
  //
  // Measurement prediction for radar
  //
  int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  
  z_pred.fill(0);
  for (int i = 0; i < n_sig_; i++) {
    //transform sigma points into measurement space
    VectorXd sigpt = Xsig_pred_.col(i);
    float px   = sigpt(0);
    float py   = sigpt(1);
    float v    = sigpt(2);
    float psi  = sigpt(3);
    float psid = sigpt(4);
    
    // Prevent div by 0
    float eps = 0.00001;
    if (fabs(px) < eps) {
      px = eps;
    }
    if (fabs(py) < eps) {
      py = eps;
    }
    
    float rho = sqrt(pow(px, 2) + pow(py, 2));
    float phi = atan2(py, px);
    float rhod = v * (px * cos(psi) + py * sin(psi)) / rho;
    
    Zsig.col(i) << rho, phi, rhod;
    
    // Add to mean
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate measurement covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R << pow(std_radr_, 2), 0,                   0,
       0,                 pow(std_radphi_, 2), 0,
       0,                 0,                   pow(std_radrd_, 2);
  S.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd diff = Zsig.col(i) - z_pred;
    S += weights_(i) * diff * diff.transpose();
  }
  S += R;
  
  //
  // Update from measurement:
  //
  MatrixXd T = MatrixXd(n_x_, n_z);
  
  T.fill(0);
  for (int i = 0; i < n_sig_; i++) {
    // Calculate cross correlation matrix.
    MatrixXd xdiff = Xsig_pred_.col(i) - x_;
    MatrixXd zdiff = Zsig.col(i) - z_pred;
    T += weights_(i) * xdiff * zdiff.transpose();
  }
  
  // Calculate Kalman gain K.
  MatrixXd K = T * S.inverse();
  
  // Update state mean and covariance matrix.
  MatrixXd z = meas_package.raw_measurements_;

  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
}
