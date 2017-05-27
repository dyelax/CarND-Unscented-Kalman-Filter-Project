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
  P_ = MatrixXd(n_x_, n_x_)::Identity();

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
  Xsig_pred_ = MatrixXd(n_aug_, n_sig_);
  
  // Init weights
  weights_ = VectorXd(n_aug_);
  for (int i = 0; i < n_sig_; i++) {
    if (i == 0) {
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    } else {
      weights_(i) = 1 / (2 * (lambda_ + n_aug_));
    }
  }
}

UKF::~UKF() {}

VectorXd PolarToCartesian(const VectorXd &polar_measurements) {
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
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  long long timestamp = measurement_pack.timestamp_;
  
  if (!is_initialized_) {
    previous_timestamp_ = timestamp;
    
    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      VectorXd cartesian_measurements =
      PolarToCartesian(measurement_pack.raw_measurements_);
      
      float px = cartesian_measurements(0);
      float py = cartesian_measurements(1);
      float v  = cartesian_measurements(2);
      
      x_ << px, py, v, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      x_ << measurement_pack.raw_measurements_(0),
            measurement_pack.raw_measurements_(1),
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

  // Update the process noise covariance matrix (Q).
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(measurement_pack.raw_measurements_);
  } else {
    UpdateLidar(measurement_pack.raw_measurements_);
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
  // Compute motion and noise matrices.
  px   = x_(0);
  py   = x_(1);
  v    = x_(2);
  yaw  = x_(3);
  yawd = x_(4);
  std_a2 = pow(std_a_, 2);
  std_yawdd2 = pow(std_yawdd_, 2);
  
  VectorXd noise = VectorXd(n_x_);
  float half_dt2 = 0.5 * pow(delta_t, 2);
  noise << half_dt2 * cos(yaw) * std_a2,
           half_dt2 * sin(yaw) * std_a2,
           delta_t * std_a2,
           half_dt2 * std_yawdd2,
           delta_t * std_yawdd2;
  
  VectorXd motion = VectorXd(n_x_);
  if (yawd == 0) {
    motion << v * cos(yaw) * delta_t,
              v * sin(yaw) * delta_t,
              0,
              yaw * delta_t,
              0;
  } else {
    motion << (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw)),
              (v / yawd) * (-cos(yaw + yawd * delta_t) + cos(yaw)),
              0,
              yaw * delta_t,
              0;
  }
  
  // Find sigma points.
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(6) = pow(std_a_, 2)
  
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a2, 0,
       0,      std_yawdd2;
  
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  
  
  MatrixXd Xsig = MatrixXd(n_aug_, n_sig_);
  
  
  
  
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
}
