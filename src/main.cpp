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
#include "spline.h"

using namespace Eigen;
using namespace std;
using namespace tk;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{
	auto found_null = s.find("null");
	auto b1 = s.find_first_of("[");
	auto b2 = s.find_first_of("}");
	if (found_null != string::npos)
	{
		return "";
	}
	else if (b1 != string::npos && b2 != string::npos)
	{
		return s.substr(b1, b2 - b1 + 2);
	}
	return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for (int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x, y, map_x, map_y);
		if (dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}

	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y - y), (map_x - x));

	double angle = fabs(theta - heading);
	angle = min(2 * pi() - angle, angle);

	if (angle > pi() / 4)
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
	int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp - 1;
	if (next_wp == 0)
	{
		prev_wp = maps_x.size() - 1;
	}

	double n_x = maps_x[next_wp] - maps_x[prev_wp];
	double n_y = maps_y[next_wp] - maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
	double proj_x = proj_norm * n_x;
	double proj_y = proj_norm * n_y;

	double frenet_d = distance(x_x, x_y, proj_x, proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x, center_y, x_x, x_y);
	double centerToRef = distance(center_x, center_y, proj_x, proj_y);

	if (centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for (int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
	}

	frenet_s += distance(0, 0, proj_x, proj_y);

	return {frenet_s, frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp + 1) % maps_x.size();

	double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s - maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
	double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

	double perp_heading = heading - pi() / 2;

	double x = seg_x + d * cos(perp_heading);
	double y = seg_y + d * sin(perp_heading);

	return {x, y};
}

double mhp_to_mps(double mhp)
{
	return mhp / 2.236936;
}

double mps_to_mhp(double mps)
{
	return mps * 2.236936;
}

vector<double> JMT(vector<double> start, vector<double> end, double T)
{
	MatrixXd A = MatrixXd(3, 3);
	A << T * T * T, T * T * T * T, T * T * T * T * T,
		3 * T * T, 4 * T * T * T, 5 * T * T * T * T,
		6 * T, 12 * T * T, 20 * T * T * T;

	MatrixXd B = MatrixXd(3, 1);
	B << end[0] - (start[0] + start[1] * T + .5 * start[2] * T * T),
		end[1] - (start[1] + start[2] * T),
		end[2] - start[2];

	MatrixXd Ai = A.inverse();

	MatrixXd C = Ai * B;

	vector<double> result = {start[0], start[1], .5 * start[2]};
	for (int i = 0; i < C.size(); i++)
	{
		result.push_back(C.data()[i]);
	}

	return result;
}

vector<double> calc_poly(vector<double> coeffs, double t)
{
	double s = coeffs[0] + coeffs[1] * t + coeffs[2] * t * t +
			   coeffs[3] * t * t * t + coeffs[4] * t * t * t * t + coeffs[5] * t * t * t * t * t;

	double s_d = coeffs[1] + 2 * coeffs[2] * t +
				 3 * coeffs[3] * t * t + 4 * coeffs[4] * t * t * t + 5 * coeffs[5] * t * t * t * t;

	double s_d_d = 2 * coeffs[2] + 3 * 2 * coeffs[3] * t + 4 * 3 * coeffs[4] * t * t + 5 * 4 * coeffs[5] * t * t * t;

	return {s, s_d, s_d_d};
}

void ConvertToWorld(double vx, double vy, vector<double> &x, vector<double> &y, double psi)
{
	double cos_psi = cos(psi);
	double sin_psi = sin(psi);

	for (unsigned int i = 0; i < x.size(); ++i)
	{
		double xn = x[i];
		double yn = y[i];

		x[i] = xn * cos_psi - yn * sin_psi;
		y[i] = xn * sin_psi + yn * cos_psi;

		x[i] = x[i] + vx;
		y[i] = y[i] + vy;
	}
}

void ConvertToVehicle(double vx, double vy, vector<double> &x, vector<double> &y, double psi)
{
	double cos_psi = cos(psi);
	double sin_psi = sin(psi);

	for (unsigned int i = 0; i < x.size(); ++i)
	{
		double xn = x[i] - vx;
		double yn = y[i] - vy;
		x[i] = xn * cos_psi + yn * sin_psi;
		y[i] = yn * cos_psi - xn * sin_psi;
	}
}

// state codes
#define KL 0
#define LCL 1
#define RCL 2


// constants
#define MAX_A 8

// structure to define full car state required to make decisions
struct CarState
{
	int state_code;
	int target_lane;
	int source_lane;
	double d_to_target;
	vector<double> lane_speeds;

	double s;
	
	double planned_x;
	double planned_y;
	double planned_angle;
	double planned_speed;
	

	vector<double> points_x;
	vector<double> points_y;

	vector<double> planned_points_x;
	vector<double> planned_points_y;
};

struct Car
{
	double x;
	double y;
	double s;
	double d;
	int lane;
	double vx;
	double vy;
	double speed;
	double planned_s;
};

struct SensorFusion
{
	vector<Car> forward_cars;
	vector<Car> back_cars;
	vector<double> lane_speeds;
};


// transform sensor fusion data - c
SensorFusion get_sensor_fusion(nlohmann::json sensor_fusion, double car_s, double car_s_planned, int path_size,
							   double max_speed, double car_detect_range, double speed_detect_range)
{
	vector<double> lane_speeds = {max_speed, max_speed, max_speed};
	vector<double> lane_closest_cars_ds = {10000, 10000, 10000};

	SensorFusion res;

	for (unsigned int i = 0; i < sensor_fusion.size(); i++)
	{
		double d = sensor_fusion[i][6];
		double s = sensor_fusion[i][5];
		double x = sensor_fusion[i][1];
		double y = sensor_fusion[i][2];
		double vx = sensor_fusion[i][3];
		double vy = sensor_fusion[i][4];

		double speed = sqrt(vx * vx + vy * vy);

		int lane = (int)floor(d / 4.0);
		double ds = s - car_s;
		double s_planned = s + path_size * 0.02 * speed;
		double ds_planned = s_planned - car_s_planned;

		if (ds_planned < speed_detect_range && ds_planned > 0)
		{
			if (ds_planned > 0 && lane >= 0 and lane < 3)
			{
				if (ds_planned < lane_closest_cars_ds[lane])
				{
					if (ds_planned < speed_detect_range - speed_detect_range * 0.2)
					{
						lane_speeds[lane] = speed - speed * 0.2;
					}
					else
					{
						lane_speeds[lane] = speed;
					}

					lane_closest_cars_ds[lane] = ds_planned;
				}
			}
		}

		if (fabs(ds_planned) < car_detect_range)
		{
			Car c = {
				x,
				y,
				s,
				d,
				lane,
				vx,
				vy,
				speed,
				s_planned};

			if (s_planned > car_s_planned)
			{
				res.forward_cars.push_back(c);
			}
			else
			{
				res.back_cars.push_back(c);
			}
		};
	}

	res.lane_speeds = lane_speeds;

	return res;
}

vector<int> successor_states(CarState current_state, map<int, vector<int>> transitions)
{
	vector<int> res;

	vector<int> states = transitions[current_state.state_code];
	
	for (int i=0; i < states.size(); i++)
	{
		int state = states[i];

		// check if car is out of road - in the real world it might go to cost function
		// as it might be cheeper than hit another car - but for the project for simplicity we filter
		// actions which leads out of the road.
		if (current_state.source_lane == 0 && state == LCL)
		{
			continue;
		}

		if (current_state.source_lane == 2 && state == RCL)
		{
			continue;
		}

		// prevent to switch states in the middle of lane changing, idealy should go to cost functions as well
		// as sometimes it is valid behaviour due to safety
		if (current_state.state_code != KL && state == KL && current_state.source_lane != current_state.target_lane)
		{
			continue;
		}

		res.push_back(state);
	}

	return res;
}

// penalty fpr changing lane  - prevents unnessesary lane changing behaviour
double line_change_cost(const CarState &current_state, const CarState &target_state, const SensorFusion &sensor_fusion,
	const vector<double> & map_waypoints_x,
	const vector<double> & map_waypoints_y,
	const vector<double> & map_waypoints_s)
{
	if (current_state.target_lane != target_state.target_lane)
	{
		double total_cost = 0.009; //constant penalty to prevent changing lines in the empty road.
		
		// dinamic penalty to prevent changing line to the lne wehre car is closer, even if it is faster -
		// AI cars have weried acceleration behaviour 

		double min_source = 9999999999999;
		double min_target = 9999999999999;


		for (int j = 0; j < sensor_fusion.forward_cars.size(); j++)
		{
			Car car = sensor_fusion.forward_cars[j];
			double car_s = car.planned_s;
			double car_d = car.d;

			if (car.lane == current_state.target_lane)
			{
				if (min_source > car_s)
				{
					min_source = car_s;
				}
			}
			else if (car.lane == target_state.target_lane)
			{
				if (min_target > car_s)
				{
					min_target = car_s;
				}
			}
		}

		if (min_target < min_source)
		{
			total_cost += 0.1;
		}

		return total_cost;
	}
	else
	{
		return 0;
	}
}

// penalty for selecting line with low speed
double speed_cost(const CarState &current_state, const CarState &target_state, const SensorFusion &sensor_fusion,
	const vector<double> & map_waypoints_x,
	const vector<double> & map_waypoints_y,
	const vector<double> & map_waypoints_s)
{
	return 1.0 / target_state.lane_speeds[target_state.target_lane];
}

// huge penalty for unsafe driving - just check futher trajectory poitns and predict how far cars will be enar these points
// if car to close - assign penalty
double safety_cost(const CarState &current_state, const CarState &target_state, const SensorFusion &sensor_fusion,
	const vector<double> & map_waypoints_x,
	const vector<double> & map_waypoints_y,
	const vector<double> & map_waypoints_s)
{

	double distance_threshold = 40.0;

	// calculate penalty for being close to the car in front.
	for (int i = 0; i < target_state.planned_points_x.size(); i += 3)
	{
		double x = target_state.planned_points_x[i];
		double y = target_state.planned_points_y[i];

		for (int j = 0; j < sensor_fusion.forward_cars.size(); j++)
		{
			Car car = sensor_fusion.forward_cars[j];
			double car_s = car.planned_s + i * 0.02 * car.speed * 0.8; // add risk that car is breaking 
			double car_d = car.d;

			vector<double> xy = getXY(car_s, car_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
			double car_x = xy[0]; 
			double car_y = xy[1];

			double dist = (car_x - x) * (car_x - x) + (car_y - y) * (car_y - y);

			
			if (dist < distance_threshold && (car.lane == target_state.target_lane || car.lane == current_state.target_lane))
			{
				return 1000.0 * 1.0 / dist;
			} 
		}

		// calculate penalty for being close to the car in back.
		for (int j = 0; j < sensor_fusion.back_cars.size(); j++)
		{
			Car car = sensor_fusion.back_cars[j];
			double car_s = car.planned_s + i * 0.02 * car.speed * 2; // multiple 2 to add risk that car will be accelerating 
			double car_d = car.d;

			vector<double> xy = getXY(car_s, car_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
			double car_x = xy[0]; 
			double car_y = xy[1];

			double dist = (car_x - x) * (car_x - x) + (car_y - y) * (car_y - y);

			//cout << "Distance to car: " << dist << endl; 

			if (dist < distance_threshold * 2 && (car.lane == target_state.target_lane))
			{
				return 1000.0 * 1.0 / dist;
			} 
		}
	}

	return 0.0;
}

// generate trajectory which transfer car from current state to target state
CarState generate_trajectory(const CarState &current_state, int state_code, 
	vector<double> & map_waypoints_x,
	vector<double> & map_waypoints_y,
	vector<double> & map_waypoints_s)
{
	CarState res;
	vector<double>  points_x = current_state.points_x;
	vector<double>  points_y = current_state.points_y;
	
	// extract car position at speed the end of the previous trajectory
	double pos_x = current_state.planned_x;
	double pos_y = current_state.planned_y;
	double angle = current_state.planned_angle;
	
	double speed_end_path = current_state.planned_speed;

	vector<double> planned_points_x;
	vector<double> planned_points_y;
	
	// extract source lane and target lane according to state transition
	res.target_lane = current_state.source_lane;
	res.source_lane = current_state.source_lane;

	double prev_size = 0;

	if (state_code == LCL)
	{
		res.target_lane -= 1;
	}

	if (state_code == RCL)
	{
		res.target_lane += 1;
	}

	int car_lane = res.target_lane;

	// get target speed of the desired lane
	double target_speed = current_state.lane_speeds[car_lane];


	// build spline which predict smooth trajectory along the road to the target lane.
	double target_d = 2.0 + 4.0 * car_lane;
	double s_scale = 60;

	double car_s = current_state.s;

	// prepare poitns for hte spline
	for (int i = 1; i < 4; i++)
	{
		double new_s = car_s + s_scale * i;
		double new_d = target_d;
		vector<double> xy = getXY(new_s, new_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

		points_x.push_back(xy[0]);
		points_y.push_back(xy[1]);
	}

	// convert points to vehilce coordinates 
	ConvertToVehicle(pos_x, pos_y, points_x, points_y, angle);

	// build spline along the points
	spline s;
	s.set_points(points_x, points_y);


	double target_x = 30.0;
	double target_y = s(target_x);
	double dist = sqrt(target_x * target_x + target_y * target_y);

	double x_start = 0;

	vector<double> points_x_v;
	vector<double> points_y_v;
	
	double next_speed = target_speed;
	
	// check if we need to brake to reach target speed
	bool isBrake = (target_speed - speed_end_path) < 0;

	for (int i = 0; i < 50; i++)
	{
		// check if acceleration exeed maximum alowed - if so adjust target speed to the value 
		// which is possible to reach with maximum acceleration 
		double acc = fabs(target_speed - speed_end_path) / 0.02;

		if (acc > MAX_A)
		{
			if (isBrake)
			{
				next_speed = speed_end_path - MAX_A * 0.02;
			}
			else
			{
				next_speed = speed_end_path + MAX_A * 0.02;
			}

			speed_end_path = next_speed;
		}

		// determine step of the point - so we can reach next planend speed
		double N = dist / (.02 * next_speed);

		// calculate next desired trajectory point
		double x = x_start + target_x / N;
		x_start = x;
		double y = s(x);

		// add planned points to the trajectory
		points_x_v.push_back(x);
		points_y_v.push_back(y);
	}

	// remember we are in car coordinates so need go back to world coordinates
	ConvertToWorld(pos_x, pos_y, points_x_v, points_y_v, angle);

	// update target state with trajectory points
	for (unsigned int i = 0; i < points_x_v.size(); ++i)
	{
		planned_points_x.push_back(points_x_v[i]);
		planned_points_y.push_back(points_y_v[i]);
	}

	// update target state with data from current state
	res.planned_points_x = planned_points_x;
	res.planned_points_y = planned_points_y;
	res.lane_speeds = current_state.lane_speeds;
	res.planned_angle = current_state.planned_angle;
	res.planned_x = current_state.planned_x;
	res.planned_y = current_state.planned_y;
	res.s = current_state.s;
	res.state_code = state_code;
	
	return res;
}

CarState transition_function(CarState current_state, map<int, vector<int>> transitions, SensorFusion sensor_fusion,
	vector<double> & map_waypoints_x,
	vector<double> & map_waypoints_y,
	vector<double> & map_waypoints_s)
{
	// this is canonical transition function. It checks all possible transition states,
	// calculates cost functions 
	// selecat target state with smallest cast
	//	# only consider states which can be reached from current state.
	vector<int> possible_successor_states = successor_states(current_state, transitions);

	//	# keep track of the total cost of each state.
	map<int, double> costs;
	map<int, CarState> states;

	vector<function<double(const CarState &, const CarState &, const SensorFusion &, const vector<double> &, const vector<double> &, const vector<double> &)>> cf_list = {line_change_cost, speed_cost, safety_cost};
	vector<double> weight_list = {1, 1, 1};

	//cout << possible_successor_states.size() << "    " << endl;

	for (int i = 0; i < possible_successor_states.size(); i++)
	{
		int state_code = possible_successor_states[i];

		// # generate trajectory
		CarState target_state = generate_trajectory(current_state, state_code,
						map_waypoints_x,
						map_waypoints_y,
						map_waypoints_s);

		double cost_for_state = 0;

		// calculate cost functions for the trajectory
		for (int j = 0; j < cf_list.size(); j++)
		{
			float new_cost = weight_list[j] * cf_list[j](current_state, target_state, sensor_fusion, map_waypoints_x, map_waypoints_y, map_waypoints_s);
			cost_for_state += new_cost;
		}

		// trace desision
		cout << "considering target state - " << state_code << endl;
		cout << "current lane - " << current_state.target_lane << ", target lane - " << target_state.target_lane << endl;
		cout << "cost for the target state - " << cost_for_state << endl;

		costs[state_code] = cost_for_state;
		states[state_code] = target_state;
	}


	// select target trajectory with smallest cost
	double min_cost = 999999999;
	CarState min_state;

	map<int, double>::iterator it;
	for (it = costs.begin(); it != costs.end(); it++)
	{
		double cost = it->second;
		double code = it->first;

		if (cost < min_cost)
		{
			min_cost = cost;
			min_state = states[code];
		}
	}

	// trace desision
	cout << "min_state:  lane - " << min_state.target_lane << "   code - " << min_state.state_code << endl;  

	return min_state;
}


int main()
{
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
	while (getline(in_map_, line))
	{
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

	//posible transitions
	map<int, vector<int>> transitions;
	transitions[KL] = {KL, RCL, LCL};
	transitions[RCL] = {KL, RCL};
	transitions[LCL] = {KL, LCL};

	// set current car state initial information
	CarState current_state = {
		KL,
		1,
		1,
		1,
		{0, 0, 0}};

	//clock_t time_start = clock();

	h.onMessage([&current_state, &transitions, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
																																						uWS::OpCode opCode) {
		// "42" at the start of the message means there's a websocket message event.
		// The 4 signifies a websocket message
		// The 2 signifies a websocket event
		//auto sdata = string(data).substr(0, length);
		//cout << sdata << endl;
		if (length && length > 2 && data[0] == '4' && data[1] == '2')
		{

			auto s = hasData(data);

			if (s != "")
			{
				auto j = json::parse(s);

				string event = j[0].get<string>();

				if (event == "telemetry")
				{
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

					//setup constants
					double max_speed = mhp_to_mps(49.0);
			
					vector<double> next_x_vals;
					vector<double> next_y_vals;

					// set the desired target speed
					double car_mps = mhp_to_mps(car_speed);
					double angle = deg2rad(car_yaw);

					double current_a = 0;

					int path_size = previous_path_x.size();

					double pos_x;
					double pos_y;
					double speed_x;
					double speed_y;

					double ax;
					double ay;

					vector<double> points_x;
					vector<double> points_y;
					double end_s;
					double end_d;

					// we are going to use previous calculated points as 
					// discussed in walkthrough - this guarantee smooth path created
					double prev_size = previous_path_x.size();
					double speed_end_path = car_mps;

					if (prev_size >= 2)
					{
						double x_1 = previous_path_x[prev_size - 1];
						double y_1 = previous_path_y[prev_size - 1];

						double x_2 = previous_path_x[prev_size - 2];
						double y_2 = previous_path_y[prev_size - 2];

						pos_x = x_1;
						pos_y = y_1;

						angle = atan2(pos_y - y_2, pos_x - x_2);

						points_x.push_back(x_2);
						points_y.push_back(y_2);

						points_x.push_back(pos_x);
						points_y.push_back(pos_y);

						speed_end_path = sqrt((x_1 - x_2) * (x_1 - x_2) + (y_1 - y_2) * (y_1 - y_2)) / 0.02;

						for (int i = 0; i < path_size; i++)
						{
							next_x_vals.push_back(previous_path_x[i]);
							next_y_vals.push_back(previous_path_y[i]);
						}
					}
					else
					{
						pos_x = car_x;
						pos_y = car_y;

						double x_2 = pos_x - cos(car_yaw);
						double y_2 = pos_y - sin(car_yaw);

						points_x.push_back(x_2);
						points_y.push_back(y_2);

						points_x.push_back(pos_x);
						points_y.push_back(pos_y);
					}


					// sensor fusion section
					// prepare sensor fusion data in a way it can be used for path generation
					// determine closest cars and lane speeds
					vector<double> sd = getFrenet(pos_x, pos_y, angle, map_waypoints_x, map_waypoints_y);
					double car_s_planned = sd[0];
					double car_d_planned = sd[1];

					SensorFusion sensor_fusion_data = get_sensor_fusion(sensor_fusion, car_s, car_s_planned, prev_size,
																		max_speed, 80.0, 30.0);
					
					// collect all nessesary informatiion to the current car state 
					// including sensor fusion data, car position and planned car position at the end ot the remaining
					// trajectory points
					current_state.lane_speeds = sensor_fusion_data.lane_speeds;
					current_state.s = car_s;
					current_state.points_x = points_x;
					current_state.points_y = points_y;

					current_state.planned_x = pos_x;
					current_state.planned_y = pos_y;
					current_state.planned_angle = angle;
					current_state.planned_speed = speed_end_path;

					// how close car to the center of taget lane - used to determine if next maneuver is allowed
					current_state.d_to_target = fabs(2 + current_state.target_lane * 4 - car_d_planned);

					// destination reached to target lane - so swiching lane in state
					if (current_state.d_to_target < 0.5) {
						current_state.source_lane = current_state.target_lane;
					}

					// trace that current cycle is completed
					cout << endl;
					cout << "Next cycle: current_lane - " <<  (car_d_planned - 2.0) / 4.0 << " state - " << current_state.state_code <<  endl;
					cout << endl;
					
					// determine target state based on the avalaible trajectories and cost functions
					CarState target_state = transition_function(current_state, transitions, sensor_fusion_data,
						map_waypoints_x,
						map_waypoints_y,
						map_waypoints_s);

					// update current state and points according to the selected target state
					current_state.state_code = target_state.state_code;
					current_state.target_lane = target_state.target_lane;

					for (unsigned int i = 0; i < 50 - prev_size; ++i)
					{
						next_x_vals.push_back(target_state.planned_points_x[i]);
						next_y_vals.push_back(target_state.planned_points_y[i]);
					}

				
					msgJson["next_x"] = next_x_vals;
					msgJson["next_y"] = next_y_vals;

					auto msg = "42[\"control\"," + msgJson.dump() + "]";

					//this_thread::sleep_for(chrono::milliseconds(1000));
					ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
				}
			}
			else
			{
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
		if (req.getUrl().valueLength == 1)
		{
			res->end(s.data(), s.length());
		}
		else
		{
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
	if (h.listen(port))
	{
		std::cout << "Listening to port " << port << std::endl;
	}
	else
	{
		std::cerr << "Failed to listen to port" << std::endl;
		return -1;
	}
	h.run();
}
