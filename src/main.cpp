#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"

#include "json.hpp"


using namespace Eigen;
using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

double mhp_to_mps(double mhp)
{
	return mhp / 2.236936;
}

double mps_to_mhp(double mps)
{
	return mps * 2.236936;
}

vector<double> JMT(vector< double> start, vector <double> end, double T)
{
    /*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficent in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    */
    
    MatrixXd A = MatrixXd(3, 3);
	A << T*T*T, T*T*T*T, T*T*T*T*T,
			    3*T*T, 4*T*T*T,5*T*T*T*T,
			    6*T, 12*T*T, 20*T*T*T;
		
	MatrixXd B = MatrixXd(3,1);	    
	B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
			    end[1]-(start[1]+start[2]*T),
			    end[2]-start[2];
			    
	MatrixXd Ai = A.inverse();
	
	MatrixXd C = Ai*B;
	
	vector <double> result = {start[0], start[1], .5*start[2]};
	for(int i = 0; i < C.size(); i++)
	{
	    result.push_back(C.data()[i]);
	}
	
    return result;
    
}

vector<double> calc_poly(vector< double> coeffs, double t)
{
	double s = coeffs[0] + coeffs[1] * t + coeffs[2] * t * t +
			coeffs[3] * t * t * t + coeffs[4] * t * t * t * t + coeffs[5] * t * t * t * t * t;

	double s_d = coeffs[1] + 2 * coeffs[2] * t +
			3 * coeffs[3] * t * t  +  4 * coeffs[4] *  t * t * t + 5 * coeffs[5] * t * t * t * t;

	double s_d_d = 2 * coeffs[2]  + 3 * 2 * coeffs[3] * t  +  4 * 3 * coeffs[4]  * t * t + 5 * 4 * coeffs[5] * t * t * t;

	return { s, s_d, s_d_d };
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }


  // state

  string current_state;
  int car_lane = 1;



  h.onMessage([&car_lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

						// to do code here


          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

			// set the desired target speed
			double max_speed = mhp_to_mps(49.0);


			int line_number = 1;


			double car_mps = mhp_to_mps(car_speed);

			double delta_t = 1.5;

			//if (car_mps < 5)
			//{
			//	delta_t = 2.0;
			//}
			//else
			//{
			//	delta_t = car_mps / 25.0;
			//}

			//double meters_to_maneur = 100.0;

			//double expected_average_speed =



			double target_speed = car_mps + 9.2 * delta_t;

			if (target_speed > max_speed) 
			{
				target_speed = max_speed;
			}

			double angle = deg2rad(car_yaw);

			double current_a = 0;
			double prev_size = previous_path_x.size();

			if (previous_path_x.size() >= 2)
			{
				double x_1 = previous_path_x[prev_size - 1];
				double y_1 = previous_path_y[prev_size - 1];
				
				double x_2 = previous_path_x[prev_size - 2];
				double y_2 = previous_path_y[prev_size - 2];


				double x_n1 = previous_path_x[0];
				double y_n1 = previous_path_y[0];
				
				double x_n2 = previous_path_x[1];
				double y_n2 = previous_path_y[1];


				//cout << "car: " << car_x << "   " << car_y << endl;
				//cout << "prev1: " << x_1 << "   " << y_1 << endl;
				//cout << "prev2: " << x_2 << "   " << y_2 << endl;
				//cout << "next1: " << x_n1 << "   " << y_n1 << endl;
				//cout << "next2: " << x_n2 << "   " << y_n2 << endl;

				double ax = x_n2 - 2 * x_n1 + car_x;
				double ay = y_n2 - 2 * y_n1 + car_y;
				
				current_a = sqrt(ax * ax + ay * ay);
				cout << "a: " << current_a << endl;

			}


			if (true)
			{
				 int path_size = previous_path_x.size();
				 double pos_x;
          		 double pos_y;
         		 double angle;
			 	 double speed_x;
				 double speed_y;

				double ax;
				double ay;
				
				
				//if (path_size > 30)
				//{
				//	path_size = 30;
				//}


				for(int i = 0; i < path_size; i++)
				{
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}

				if(path_size == 0)
				{
					pos_x = car_x;
					pos_y = car_y;
					angle = deg2rad(car_yaw);
					speed_x = car_mps * cos(angle);
					speed_y = car_mps * sin(angle);
					ax = 0;
					ay = 0;
				}
				else
				{
					pos_x = previous_path_x[path_size-1];
					pos_y = previous_path_y[path_size-1];

					double pos_x2 = previous_path_x[path_size-2];
					double pos_y2 = previous_path_y[path_size-2];
					angle = atan2(pos_y-pos_y2,pos_x-pos_x2);

			     	double pos_x3 = previous_path_x[path_size-3];
					double pos_y3 = previous_path_y[path_size-3];

					ax = (pos_x3 - 2 * pos_x2 + pos_x) / 0.02 / 0.02;
					ay = (pos_y3 - 2 * pos_y2 + pos_y) / 0.02 / 0.02;

					speed_x = (pos_x - pos_x2) / 0.02;
					speed_y = (pos_y - pos_y2) / 0.02;
		
				}

				vector<double> frenet = getFrenet(pos_x, pos_y, angle, map_waypoints_x, map_waypoints_y);
				
				double current_s = frenet[0];
				double current_d = frenet[1];
				
				double target_s = current_s + 30;
				double target_d = 2 + 4 * car_lane;
				

				//double speed_x = car_mps * cos(angle);
				//double speed_y = car_mps * sin(angle);

				double target_speed_x = target_speed * cos(angle);
				double target_speed_y = target_speed * sin(angle);
				 
				cout << "target speed x: " << target_speed_x << "target speed y: " << target_speed_y << endl;
				cout << "current speed x: " << speed_x << "current speed y: " << speed_y << endl;
				


				vector<double> car_xy = { pos_x, pos_y };
				vector<double> target_xy = getXY(target_s, target_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

				vector<double> coefs_x = JMT({car_xy[0], speed_x, ax}, {target_xy[0], target_speed_x, 0}, delta_t);
				vector<double> coefs_y = JMT({car_xy[1], speed_y, ay}, {target_xy[1], target_speed_y, 0}, delta_t);
				

				for(int i = 0; i < 50-path_size; i++)
				{
					vector<double> x_ = calc_poly(coefs_x, 0.02 * (i + 1));
					vector<double> y_ = calc_poly(coefs_y, 0.02 * (i + 1));

					cout << "point speed x: " << x_[1] << "point speed y: " << y_[1] << endl;
					
					//cout << "S: " << new_s << "  speed: " << s_s[1] << endl;
					//cout << "curr: " << xy[0] << " " << xy[1] << endl;
				
					next_x_vals.push_back(x_[0]);
					next_y_vals.push_back(y_[0]);
				}


				// state
				double state_x = 0;
				double state_x_d = 0;
				double state_x_dd = 0;
				
				double state_y = 0;
				double state_y_d = 0;
				double state_y_dd = 0;
				

				double state_x_f = (car_mps  + target_speed) / 2.0  * delta_t;
				double state_x_d_f = target_speed;
				double state_y_dd_f = 0;
				
				double state_d_f = 0;
				double state_d_d_f = 0;
				double state_d_dd_f = 0;
				

				
			
				}

			else
			{




				// state
				double state_s = 0;
				double state_s_d = car_mps;
				double state_s_dd = current_a;
				
				double state_d = 0;
				double state_d_d = 0;
				double state_d_dd = 0;
				

				double state_s_f = (car_mps  + target_speed) / 2.0  * delta_t;
				double state_s_d_f = target_speed;
				double state_s_dd_f = 0;
				
				double state_d_f = 0;
				double state_d_d_f = 0;
				double state_d_dd_f = 0;
				

				vector<double> coefs_s = JMT({state_s, state_s_d, state_s_dd}, {state_s_f, state_s_d_f, state_s_dd_f}, delta_t);
				vector<double> coefs_d = JMT({state_d, state_d_d, state_d_dd}, {state_d_f, state_d_d_f, state_d_dd_f}, delta_t);


				cout << "Speeds:" <<  car_mps << "           " << target_speed << endl;
				cout << "S:" <<  car_s << "           " << car_s + state_s_f << endl;

				double dist_inc = 1.0 * target_speed / 50.0;

				for(int i = 0; i < 500; i++)
				{
					vector<double> s_s = calc_poly(coefs_s, 0.02 * (i + 1));
					double new_s = car_s + s_s[0];
					double new_d = car_d;
					vector<double> xy = getXY(new_s, new_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

					//cout << "S: " << new_s << "  speed: " << s_s[1] << endl;
					//cout << "curr: " << xy[0] << " " << xy[1] << endl;
				
					next_x_vals.push_back(xy[0]);
					next_y_vals.push_back(xy[1]);
				}

			}
          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	
						//double dist_inc = 1.0 * target_speed / 50.0;
						//for(int i = 0; i < 50; i++)
						//{
					//				double new_s = car_s + dist_inc*i;
				//					double new_d = car_d;
			//						vector<double> xy = getXY(new_s, new_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

//									next_x_vals.push_back(xy[0]);
//									next_y_vals.push_back(xy[1]);
//						}

						//double dist_inc = 0.1;
						//for(int i = 0; i < 50; i++)
						//{
						//			next_x_vals.push_back(car_x+(dist_inc*i)*cos(deg2rad(car_yaw)));
						//			next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad(car_yaw)));
						//}
						
						//	double pos_x;
						//	double pos_y;
						//	double angle;
						//	int path_size = previous_path_x.size();

						//	for(int i = 0; i < path_size; i++)
						//	{
						//			next_x_vals.push_back(previous_path_x[i]);
						//			next_y_vals.push_back(previous_path_y[i]);
						//	}

						//	if(path_size == 0)
						//	{
						//			pos_x = car_x;
						//			pos_y = car_y;
						//			angle = deg2rad(car_yaw);
						//	}
						//	else
						//	{
						//			pos_x = previous_path_x[path_size-1];
						//			pos_y = previous_path_y[path_size-1];

						//			double pos_x2 = previous_path_x[path_size-2];
						//			double pos_y2 = previous_path_y[path_size-2];
						//			angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
						//	}

						//	double dist_inc = 0.5;
						//	for(int i = 0; i < 50-path_size; i++)
						//	{    
						//			next_x_vals.push_back(pos_x+(dist_inc)*cos(angle+(i+1)*(pi()/100)));
						//			next_y_vals.push_back(pos_y+(dist_inc)*sin(angle+(i+1)*(pi()/100)));
						//			pos_x += (dist_inc)*cos(angle+(i+1)*(pi()/100));
						//			pos_y += (dist_inc)*sin(angle+(i+1)*(pi()/100));
						//	}


          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;











						
						//
						//
						//
						//
						
						msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

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


void transition_function(string current_state)
{
	//def transition_function(predictions, current_fsm_state, current_pose, cost_functions, weights):
	//		# only consider states which can be reached from current FSM state.
	//		possible_successor_states = successor_states(current_fsm_state)

	//			# keep track of the total cost of each state.
	//		costs = []
	//		for state in possible_successor_states:
	//				# generate a rough idea of what trajectory we would
	//				# follow IF we chose this state.
	//				trajectory_for_state = generate_trajectory(state, current_pose, predictions)
	//
	//				# calculate the "cost" associated with that trajectory.
	//				cost_for_state = 0
	//				for i in range(len(cost_functions)) :
	//						# apply each cost function to the generated trajectory
	//						cost_function = cost_functions[i]
	//						cost_for_cost_function = cost_function(trajectory_for_state, predictions)
	//
	//						# multiply the cost by the associated weight
	//						weight = weights[i]
	//						cost_for_state += weight * cost_for_cost_function
	//				costs.append({'state' : state, 'cost' : cost_for_state})
	//
	//		# Find the minimum cost state.
	//		best_next_state = None
	//		min_cost = 9999999
	//		for i in range(len(possible_successor_states)):
	//				state = possible_successor_states[i]
	//				cost  = costs[i]
	//				if cost < min_cost:
	//						min_cost = cost
	//						best_next_state = state 

   // return best_next_state

}
