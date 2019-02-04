#include <math.h>
#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include "FusionEKF_fp.h"
#include "tools_fp.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::string;
using std::vector;
using ekf::Params;
using ekf::State; 

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("]");
  if (found_null != string::npos) {
    return "";
  }
  else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

int main() {
  uWS::Hub h;

  // used to compute the RMSE later  
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  const Params& params = ekf::setup_params();
  bool initialized = false;

  ekf::State state;

  h.onMessage([&params, &initialized, &state, &estimations, &ground_truth]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, 
               uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      auto s = hasData(string(data));

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          string sensor_measurement = j[1]["sensor_measurement"];
          
          MeasurementPackage meas_package;
          std::istringstream iss(sensor_measurement);
          
          long long timestamp;

          // reads first element from the current line
          string sensor_type;
          iss >> sensor_type;

          if (sensor_type.compare("L") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::LASER;
            meas_package.raw_measurements_ = VectorXd(2);
            double px, py;            
            iss >> px >> py;            
            meas_package.raw_measurements_ << px, py;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
          } else if (sensor_type.compare("R") == 0) {
            meas_package.sensor_type_ = MeasurementPackage::RADAR;
            meas_package.raw_measurements_ = VectorXd(3);
            double ro, theta, ro_dot;            
            iss >> ro >> theta >> ro_dot;
            meas_package.raw_measurements_ << ro, theta, ro_dot;
            iss >> timestamp;
            meas_package.timestamp_ = timestamp;
          }

          double x_gt, y_gt, vx_gt, vy_gt;
          iss >> x_gt >> y_gt >> vx_gt >> vy_gt;

          VectorXd gt_values(4);
          gt_values << x_gt, y_gt, vx_gt, vy_gt;
          ground_truth.push_back(gt_values);
          
          // Call ProcessMeasurement(meas_package) for Kalman filter
          if( !initialized ) {
            state = ekf::initial_state( meas_package ); 
            initialized = true; 
          } else {
            state = ekf::proc_measurement(state, meas_package, params);            
          }

          // Push the current estimated x,y positon from the Kalman filter's 
          //   state vector
          
          double p_x = state.x_(0),
                 p_y = state.x_(1),
                 v_x = state.x_(2),
                 v_y = state.x_(3);

          VectorXd estimate(4);
          estimate << p_x, p_y, v_x, v_y;
        
          estimations.push_back(estimate);

          VectorXd RMSE = tools::CalculateRMSE(estimations, ground_truth);

          // TODO: take out 
          // cout << RMSE << endl; 

          json msgJson;
          msgJson["estimate_x"] = p_x;
          msgJson["estimate_y"] = p_y;
          msgJson["rmse_x"] =  RMSE(0);
          msgJson["rmse_y"] =  RMSE(1);
          msgJson["rmse_vx"] = RMSE(2);
          msgJson["rmse_vy"] = RMSE(3);
          auto msg = "42[\"estimate_marker\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }  // end "telemetry" if

      } else {
        string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }  // end websocket message if

  }); // end h.onMessage

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, 
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  
  h.run();
}