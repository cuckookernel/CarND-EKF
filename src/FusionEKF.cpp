#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF( ) {
  
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process noises
   */
  H_laser_ << 1.0,   0, 0, 0,
                0, 1.0, 0, 0;

  MatrixXd P(4,4);
  P << 1, 0,      0,    0,
       0, 1,      0,    0,
       0, 0,     10,    0,
       0, 0,      0,   10;

  MatrixXd F = MatrixXd::Identity(4,4);

  // Process noise is set everytime a new measurement comes in, as it depends on dt 
  ekf_.InitPF( P, F );

  noise_ax = 9;
  noise_ay = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::Initialize( const MeasurementPackage &meas_pack ) {
   /**     
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    
    // ekf_.x_ << 1, 1, 1, 1;  // Mateo: Commented out
    VectorXd rm = meas_pack.raw_measurements_;

    if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {            
      float rho     = rm(0),
            phi     = rm(1), 
            rho_dot = rm(2);

      // TODO: check
      float px = rho * cos(phi), 
            py = rho * sin(phi), 
            vx = rho_dot * cos(phi),
            vy = rho_dot * sin(phi);

      ekf_.x_ = VectorXd(4);
      ekf_.x_ << px, py, vx, vy;
    }
    else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {            
      float px = rm(0), py = rm(1);
      ekf_.x_ = VectorXd(4);
      ekf_.x_ << px, py, 0, 0; 
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
} 

void FusionEKF::ProcessMeasurement(const MeasurementPackage &meas_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    Initialize( meas_pack );
    return;
  }

  /**
   * Prediction
   */
  float dt = (meas_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_pack.timestamp_; 
  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * */
  
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  /*
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
 
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4 / 4 * noise_ax, 0, dt_3/2 * noise_ax, 0,
              0, dt_4 / 4 * noise_ay, 0, dt_3/2 * noise_ay,
              dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax  , 0,
              0, dt_3 / 2 * noise_ay, 0,   dt_2 * noise_ay;

  ekf_.Predict();
  /**
   * Update
   * - Update the state and covariance matrices.
   */

  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates    
    Hj_ = tools.CalculateJacobian( ekf_.x_ );
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF( meas_pack.raw_measurements_ );    
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    // 
    ekf_.Update( meas_pack.raw_measurements_ );
  }

  // print the output
  if( verbose_ ) {
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
  }
}
