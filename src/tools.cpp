#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  long n = estimations.size();
  if( n < 1 || n != ground_truth.size() ) {
      cout << "Invalid inputs est_size=" << estimations.size() << endl;
      return rmse;
  }

  VectorXd accum(4);
  accum << 0, 0, 0, 0;
  for (int i=0; i < estimations.size(); ++i) {
    VectorXd diff = ( estimations[i] - ground_truth[i] ); 
    accum += diff.cwiseProduct( diff ); 
  }

  // take the mean and sqrt for each component
  rmse  = (accum / n).array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  double norm_p2 = px*px + py*py;
  double norm_p  = sqrt( norm_p2 );
  double norm_p3 = norm_p2 * norm_p;
  
  MatrixXd Hj(3,4);
  // check division by zero
  if( norm_p2 < 1e-6 ) {
      cout << "Division by zero" << endl;
      return Hj;
  }  
  // compute the Jacobian matrix
  double det = ( vx * py - vy * px )  / norm_p3 ; 
  Hj <<  px / norm_p , py / norm_p , 0, 0,
        -py / norm_p2, px / norm_p2, 0, 0,        
        py * det, -px * det ,  px  / norm_p, py / norm_p;
  return Hj;
}
