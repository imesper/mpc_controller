#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "helpers.h"
#include "json.hpp"
#include <chrono>
#include <iostream>
#include <math.h>
#include <string>
#include <thread>
#include <uWS/uWS.h>
#include <vector>

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

extern size_t N;
extern double dt;
extern double Lf;

int main() {
    uWS::Hub h;

    // MPC is initialized here!
    MPC mpc;

    h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                       uWS::OpCode opCode) {
        // "42" at the start of the mesdwwwwwwsage means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        string sdata = string(data).substr(0, length);
//    std::cout << sdata << std::endl;
        if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
            string s = hasData(sdata);
            if (s != "") {
                auto j = json::parse(s);
                string event = j[0].get<string>();
                if (event == "telemetry") {
                    // j[1] is the data JSON object
                    vector<double> ptsx = j[1]["ptsx"];
                    vector<double> ptsy = j[1]["ptsy"];
                    double px = j[1]["x"];
                    double py = j[1]["y"];
                    double psi = j[1]["psi"];
                    double psi_unity = j[1]["psi_unity"];
                    double v = j[1]["speed"];
                    double delta = j[1]["steering_angle"];
                    double a = j[1]["throttle"];

                    vector<double> ptsx_car;
                    vector<double> ptsy_car;

                    // Transform to Car coordinates
                    for (size_t i = 0; i < ptsx.size(); i++) {
                        double ptx_diff = ptsx[i] - px;
                        double pty_diff = ptsy[i] - py;
                        ptsx_car.push_back(ptx_diff * cos(psi) + pty_diff * sin(psi));
                        ptsy_car.push_back(pty_diff * cos(psi) - ptx_diff * sin(psi));
                    }

                    Eigen::VectorXd ptx = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
                            ptsx_car.data(), static_cast<long>(ptsx_car.size()));
                    Eigen::VectorXd pty = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(
                            ptsy_car.data(), static_cast<long>(ptsy_car.size()));

                    Eigen::VectorXd coeffs = polyfit(ptx, pty, 3);

                    /*
                     * Calculate Error - The reference frame is the car
                     */
                    // Initial state.
                    double x0 = 0;
                    double y0 = 0;
                    double psi0 = 0;
                    double cte0 = coeffs[0];
                    double epsi0 = -atan(coeffs[1]);

                    // State after delay.
                    double x_delay = x0 + ( v * cos(psi0) * dt );
                    double y_delay = y0 + ( v * sin(psi0) * dt );
                    double psi_delay = psi0 - ( v * delta * dt / Lf );
                    double v_delay = v + a * dt;
                    double cte_delay = cte0 + ( v * sin(epsi0) * dt );
                    double epsi_delay = epsi0 - ( v * atan( (3*coeffs[3]*pow(x_delay, 2) + (2*coeffs[2]*x_delay) + coeffs[1] )) * dt / Lf );


                    double psi_car = 0;
                    double psi_pred = (v/Lf) * delta * dt;
                    std::cout << "PSI: "<< epsi0  << ": " << M_PI_2 << std::endl;
                    double v_pred = v + a * dt;

                    double x = v * cos(psi_car) * dt;
                    double y = v * sin(psi_car) * dt;


                    double cte = polyeval(coeffs, x) - x;
                    //double epsi =  atan(coeffs[1]);
                    double epsi =  atan(coeffs[1] + (2*coeffs[2]*x) + (3*coeffs[3]*x*x));

                    Eigen::VectorXd state(6);

                    state << x_delay, y_delay, psi_delay, v_delay, cte_delay, epsi_delay;

                    auto vars = mpc.Solve(state, coeffs);
//          std::cout << "Steering :" << vars[6*N] << std::endl;

                    /**
                     * TODO: Calculate steering angle and throttle using MPC.
                     * Both are in between [-1, 1].
                     */
                    double steer_value;
                    double throttle_value;
                    steer_value = -vars[6*N];
                    throttle_value = vars[6*N+9];

                    json msgJson;
                    // NOTE: Remember to divide by deg2rad(25) before you send the
                    //   steering value back. Otherwise the values will be in between
                    //   [-deg2rad(25), deg2rad(25] instead of [-1, 1].
                    msgJson["steering_angle"] = steer_value;
                    msgJson["throttle"] = throttle_value;

                    // Display the MPC predicted trajectory
                    vector<double> mpc_x_vals;
                    vector<double> mpc_y_vals;

                    for (size_t i =0; i < N; i++) {
                        mpc_x_vals.push_back(vars[i]);
                        mpc_y_vals.push_back(vars[i+N]);
                    }

                    msgJson["mpc_x"] = mpc_x_vals;
                    msgJson["mpc_y"] = mpc_y_vals;

                    // Display the waypoints/reference line
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;

                    for (int i = 0; i < 250; i += 5 ) {
                        next_x_vals.push_back(i);
                        next_y_vals.push_back(polyeval(coeffs, i));
                    }

                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;

                    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
//          std::cout << msg << std::endl;
                    // Latency
                    // The purpose is to mimic real driving conditions where
                    //   the car does actuate the commands instantly.
                    //
                    // Feel free to play around with this value but should be to drive
                    //   around the track with 100ms latency.
                    //
                    // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE SUBMITTING.
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                } // end "telemetry" if
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        } // end websocket if
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
