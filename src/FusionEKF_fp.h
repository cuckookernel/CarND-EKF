#ifndef FusionEKF_FP_H_
#define FusionEKF_FP_H_

#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"
#include "tools_fp.h"
#include <limits>
#include <memory> 

namespace ekf {  

    typedef long long tstamp_t;
    typedef Eigen::Matrix<double,4,1> State_x;
    typedef Eigen::Matrix<double,4,4> Pmat;
    typedef Eigen::Matrix<double,2,4> Hlaser_t;
    typedef Eigen::Matrix<double,2,2> Rlaser_t;
    typedef Eigen::Matrix<double,3,3> Rradar_t;


    struct State {
        /* Immutable record, all members const and public */ 

        public: 
            const State_x x_;
            const Pmat P_;
            const tstamp_t timestamp_ = std::numeric_limits<long long>::min();

            State() {}; 

            State( const Eigen::VectorXd& x_ini, const Eigen::MatrixXd& P_ini, tstamp_t ts_ini )
            : x_(x_ini), P_(P_ini), timestamp_(ts_ini) {};

            State( const State& rhs ) 
            : x_(rhs.x_), P_(rhs.P_), timestamp_(rhs.timestamp_) {};

            void operator=( const State& rhs ) {
                const_cast<State_x&>(x_)  = rhs.x_;
                const_cast<Pmat&>( P_ )  = rhs.P_;
                const_cast<tstamp_t&>( timestamp_ ) = rhs.timestamp_;
            }         
            
    };
        
    
    struct Params {
        public : 
            const Hlaser_t H_laser;
            const Rlaser_t R_laser;            
            const Rradar_t R_radar;            

            const double noise_ax, noise_ay;
            bool verbose;

            Params( const Hlaser_t& H_laser_in, const Rlaser_t& R_laser_in, const Rradar_t&  R_radar_in, 
                    double noise_ax_in, double noise_ay_in ) 
                    : H_laser( H_laser_in ), R_laser( R_laser_in ), R_radar( R_radar_in ), 
                      noise_ax( noise_ax_in ), noise_ay( noise_ay_in ) {};

    };

    const Params setup_params( );

    State initial_state( const MeasurementPackage& mp );

    /**
    * Prediction Predicts the state and the state covariance
    * using the process model
    * @param delta_T Time between k and k+1 in s
    */

    State predict( const State& state, const Eigen::MatrixXd& F, const Eigen::MatrixXd& Q, tstamp_t new_ts );

    State proc_measurement( const State& state,
                               const MeasurementPackage &measurement_pack, const Params& params);

    /* implement linear KF update */ 
    State update_laser( const State& state, const Eigen::VectorXd &z, const Hlaser_t& H, const Rlaser_t& R );

    /* implement extended KF update */ 
    State update_radar( const State& state, const Eigen::VectorXd &z, const Rradar_t& R );

    void maybe_log( const Params& params, const State& state );
        
}


#endif // FusionEKF_H_
