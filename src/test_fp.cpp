
#include <iostream>
#include "FusionEKF_fp.h"
#include "tools_fp.h"
#include <chrono>

using std::vector;
using std::string;
using std::cout;  
using std::endl; 
using Eigen::VectorXd;
using Eigen::MatrixXd;
using ekf::State;
using ekf::Params;


vector<string> str_split( string a_str, char sep );
string vector_to_str( VectorXd v, int precision );

class TestState {

    public:
        void proc_line( string line  );

        State state_;        
        bool initialized_ = false; 

        vector<VectorXd> estimations_;
        vector<VectorXd> ground_truths_;

        Params params_;
        
        TestState( ) : params_( ekf::setup_params() ) {};
};

int main() {

    string input_fn = "../data/obj_pose-laser-radar-synthetic-input.txt";

    TestState st;
    st.params_.verbose = false; 
    bool compute_rmse = false; 

    // used to compute the RMSE later
        
    std::ifstream ifs( input_fn );
    string line;

    auto t0 = std::chrono::system_clock::now();

    long count = 0;
    cout<< "After opening... good: " << ifs.good() << " is_open: " << ifs.is_open() << endl;
    VectorXd RMSE;
    RMSE.setZero(); 

    while( getline(ifs, line) ) {
        // cout << "\n\n FP " << count << " : " << endl;
        st.proc_line( line );        
        if( compute_rmse ) {
            RMSE    = tools::CalculateRMSE( st.estimations_, st.ground_truths_);   
        }
        
        if (st.params_.verbose ) {
            auto n = st.estimations_.size() - 1;
            assert( st.estimations_.size() == st.ground_truths_.size() );
            VectorXd diffs = st.estimations_[n] - st.ground_truths_[n];
        
            cout << count << " :  diffs: " <<  vector_to_str(diffs, 3)             
                << "  diffs_norm: " << diffs.norm()
                << "  x,y: " << st.state_.x_(0)<< ", " << st.state_.x_(1)
                << "    RMSE: " << vector_to_str(RMSE, 3) << endl;            
        }

        count ++;
        /* if( count > 3 ) {
            break;
        }*/
    }

    auto t1 = std::chrono::system_clock::now();
    auto elapsed = (t1 - t0).count();
    cout << "elapsed: " << elapsed << endl;

    ifs.close();
    return 0;
}

void TestState::proc_line( string line  ) {

    // Create a Kalman Filter instance
    vector<string> ps = str_split( line, '\t');
    
    MeasurementPackage meas_package;
    VectorXd ground_truth = VectorXd( 4 );

    if( ps[0] == "L") { // a laser measurement 
        // sensor_type 0, x_measured 1, y_measured 2, timestamp 3,
        // x_groundtruth 4, y_groundtruth 5, vx_groundtruth 6 , vy_groundtruth 7, 
        // yaw_groundtruth, yawrate_groundtruth.
        meas_package.sensor_type_ = MeasurementPackage::LASER;
        meas_package.raw_measurements_ = VectorXd(2);
        meas_package.raw_measurements_ << stod( ps[1]), stod(ps[2]) ;        
        
        meas_package.timestamp_ = stol( ps[3] );

        ground_truth << stod(ps[4]), stod( ps[5] ), stod( ps[6] ), stod( ps[7] );

    } else if ( ps[0] == "R" ) {
        // sensor_type 0, rho_measured 1, phi_measured 2, rhodot_measured 3, timestamp 4,
        //  x_groundtruth 5, y_groundtruth 6, vx_groundtruth 7, vy_groundtruth 8, yaw_groundtruth 9, yawrate_groundtruth 10 
        meas_package.sensor_type_ = MeasurementPackage::RADAR;

        meas_package.raw_measurements_ = VectorXd(3);
        float rho = stod( ps[1] ), 
              phi = stod( ps[2] ),
              rho_dot = stod( ps[3] );
        
        meas_package.raw_measurements_ << rho, phi, rho_dot;
        meas_package.timestamp_ = stol( ps[4] ); 
        
        ground_truth << stod( ps[5] ), stod( ps[6] ), stod( ps[7] ), stod( ps[8] );

    } else {
        cout << "Unexpected line: " << line << "\nps[0]=("<<ps[0]<<")" << endl; 
        assert( false  );
    }
    ground_truths_.push_back(ground_truth);
                
    // Call ProcessMeasurement(meas_package) for Kalman filter
    if( !initialized_ ) {
        state_ = ekf::initial_state( meas_package ); 
        initialized_ = true; 
    } else {
        state_ = ekf::proc_measurement(state_, meas_package, params_);       
    }

    // Push the current estimated x,y positon from the Kalman filter's 
    //   state vector    
    estimations_.push_back( state_.x_ );
}

string vector_to_str( VectorXd v, int precision = 3 ) {
    std::stringstream ss;

    ss.precision( precision );
    ss << "[";
    int n = v.size();
    for ( auto i = 0; i< n; i++ ) {
        ss << " " << v[i] << ",";
    }
    ss << "]";
    return ss.str();
}

vector<string> str_split( string a_str, char sep ) {
    std::stringstream ss( a_str ); 
    vector<string> result;

    while( ss.good() )  {
        string substr;
        getline( ss, substr, sep );
        // printf( "substr=%s", substr.c_str() );
        result.push_back( substr );
    }

    return result;
}
