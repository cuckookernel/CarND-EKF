#include "kalman_filter.h"
#include "assert.h"
#include "math.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::InitPF( MatrixXd &P_in, MatrixXd &F_in) {  
  P_ = P_in;
  F_ = F_in;  
}


void KalmanFilter::Predict() {
  // taken verbatim from 23.14 Laser Measurements part 4
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

  // cout << "Predict:\n" << x_ << endl << "P_" << P_ << endl; 

}

void KalmanFilter::Update(const VectorXd &z) {
  // taken verbatim from 23.14 Laser Measurements part 4 
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;  
  MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  /*
  cout << "Update Linear: " << endl
       <<  "z_pred= " << z_pred.transpose() << endl
       <<  "y = " << z_pred.transpose() << endl
       <<  "H = " << H_ << endl       
       << "S =" << S << endl
       << "K =" << K << endl
       << "x_=" << x_ << endl
       << "P_ =" << P_ << endl 
       ; */
}

#include <iostream>
#include <cmath>
using std::cout; // using std::{ cout, endl }
using std::endl; 

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // z_pred = h(x)  where h is non-linear
  VectorXd z_pred = VectorXd(3);  
  
  double rho = sqrt( x_(0) * x_(0) + x_(1) * x_(1) );
  double phi = atan2( x_(1), x_(0) );

  // cout << " phi=" << phi << endl ; 
  assert ( - M_PI <= phi &&  phi <=  M_PI ); 
  // assert ( 0.0 <= phi &&  phi <= 2 * M_PI ); 
  double rho_dot = ( x_(0) * x_(2) + x_(1) * x_(3) ) / rho; 
  
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;
  y(1) = fmod( y(1) , 2 * M_PI);

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;    
  MatrixXd K = P_ * Ht * S.inverse();

  VectorXd x1 = x_;
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

  /*
  cout << "Update EKF: " << endl
      <<  "x1= " << x1.transpose() << endl
       <<  "z_pred= " << z_pred.transpose() << endl
       <<  "y = " << z_pred.transpose() << endl
       <<  "H = " << H_ << endl       
       << "S =" << S << endl
       << "K =" << K << endl
       << "x_=" << x_ << endl
       << "P_ =" << P_ << endl 
       ;*/
}
