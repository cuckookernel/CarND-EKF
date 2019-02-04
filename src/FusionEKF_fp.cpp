#include "FusionEKF_fp.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools_fp.h"

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXd;

using std::cout;
using std::endl;
using std::vector;

using ekf::State;
    
namespace ekf {

const Params setup_params( ) {
  
  // initializing matrices
  // measurement matrix laser
  Matrix<double,2,4> H_laser;   
  H_laser << 1.0,   0, 0, 0,
                0, 1.0, 0, 0;

  //measurement covariance matrix - laser
  Matrix<double,2,2> R_laser;  
  R_laser << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  Matrix<double,3,3> R_radar;  
  R_radar << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // Process noise is set everytime a new measurement comes in, as it depends on dt 
  // ekf_.InitPF( P, F );

  double noise_ax = 9;
  double noise_ay = 9;

  return Params( H_laser, R_laser, R_radar, noise_ax, noise_ay );

}

StatePtr ekf::initial_state( const MeasurementPackage &meas_pack ) {

    // TODO take P from client
    MatrixXd P(4,4);
    P << 1, 0,      0,    0,
         0, 1,      0,    0,
         0, 0,   1000,    0,
         0, 0,      0,  1000;

    // first measurement
    cout << "EKF: " << endl;
    
    // ekf_.x_ << 1, 1, 1, 1;  // Mateo: Commented out
    VectorXd rm = meas_pack.raw_measurements_;
    
    if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {            
      double rho     = rm(0),
             phi     = rm(1), 
             rho_dot = rm(2);
      
      double px = rho * cos(phi), 
             py = rho * sin(phi), 
             vx = rho_dot * cos(phi),
             vy = rho_dot * sin(phi);

      VectorXd x(4);
      x << px, py, vx, vy;
      return StatePtr( new State(x, P, meas_pack.timestamp_) ); 
    }
    else if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {            
      float px = rm(0), py = rm(1);
      VectorXd x(4);
      x << px, py, 0, 0; 

      return StatePtr( new State(x, P, meas_pack.timestamp_) ); 
    }
 
} 

StatePtr predict( const StatePtr& state, const MatrixXd& F, const MatrixXd& Q ) {
  // taken verbatim from 23.14 Laser Measurements part 4
  MatrixXd Ft = F.transpose();
  auto new_x = F * state->x_;
  auto new_P = F * state->P_ * Ft + Q;

  return StatePtr( new State(new_x, new_P, state->timestamp_) ) ;
   
}


StatePtr proc_measurement(const StatePtr& state0, const MeasurementPackage &meas_pack, const Params& params) {
  /**
   * Prediction
   */
  long long new_ts = meas_pack.timestamp_;
  float dt = ( new_ts - state0->timestamp_) / 1000000.0;
  //cout << "dt = " << dt;  
  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * */
  MatrixXd F(4,4); 
  F.setIdentity();  
  F(0, 2) = dt;
  F(1, 3) = dt;
  /*
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */ 
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  
  double n_ax = params.noise_ax, 
         n_ay = params.noise_ay; 
  Matrix<double, 4, 4 > Q;
  Q <<  dt_4 / 4 * n_ax, 0, dt_3/2 * n_ax, 0,
        0, dt_4 / 4 * n_ay, 0, dt_3/2 * n_ay,
        dt_3 / 2 * n_ax, 0, dt_2 * n_ax  , 0,
        0, dt_3 / 2 * n_ay, 0,   dt_2 * n_ay;

  auto state1 = predict( state0, F, Q );
  /**
   * Update
   * - Update the state and covariance matrices.
   */
  
  if (meas_pack.sensor_type_ == MeasurementPackage::LASER) {
    // TODO: Radar updates    
    auto state2 = update_laser( state1, meas_pack.raw_measurements_, params.H_laser, params.R_laser );
    maybe_log( params, state2 );    
    return state2;    
  } else {
    // 
    auto state2 = update_radar( state1, meas_pack.raw_measurements_, params.R_radar );    
    maybe_log( params, state2 );    
    return state2;
  }
  
}

void maybe_log( const Params& params, const StatePtr& state ) {
  // print the output
    if( params.verbose ) {
      cout << "x_ = " << state->x_ << endl;
      cout << "P_ = " << state->P_ << endl;
    }
}

StatePtr update_laser(const StatePtr& state, const VectorXd &z, 
                      const Hlaser_t& H, const Rlaser_t& R, tstamp_t new_ts) {
                        
  // taken verbatim from 23.14 Laser Measurements part 4 
  const VectorXd& x1 = state->x_; 

  VectorXd z_pred = H * x1;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * state->P_ * Ht + R;  
  MatrixXd K = state->P_ * Ht * S.inverse();

  //new estimate  
  // long x_size = x1.size();
  Matrix<double, 4, 4> I;
  I.setIdentity(); //  = MatrixXd::Identity(x_size, x_size);

  auto new_x = x1 + (K * y);
  auto new_P = (I - K * H) * state->P_;
  return StatePtr( new State( new_x, new_P, new_ts ) ); 
                
}

StatePtr update_radar(const StatePtr& state, const VectorXd &z, const Rradar_t& R, tstamp_t new_ts ) {
  // z_pred = h(x)  where h is non-linear
  const VectorXd & x1 = state->x_;
  const Matrix<double,4,4> P = state->P_;

  MatrixXd H = tools::CalculateJacobian( x1 );
  
  double rho = sqrt( x1(0) * x1(0) + x1(1) * x1(1) );
  double phi = atan2( x1(1), x1(0) );
  
  // cout << " phi=" << phi << endl ; 
  assert ( - M_PI <= phi &&  phi <=  M_PI ); 
  // assert ( 0.0 <= phi &&  phi <= 2 * M_PI ); 
  double rho_dot = ( x1(0) * x1(2) + x1(1) * x1(3) ) / rho; 

  VectorXd z_pred = VectorXd(3);  
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;
  y(1) = fmod( y(1) , 2 * M_PI);

  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P * Ht + R;    
  MatrixXd K = P * Ht * S.inverse();

  Matrix<double, 4, 4> I;
  I.setIdentity(); //  = MatrixXd::Identity(x_size, x_size);

  auto new_x =  x1 + (K * y);
  auto new_P = (I - K * H) * state->P_;
  
  return StatePtr( new State( new_x, new_P, new_ts ) );   
}

}  

