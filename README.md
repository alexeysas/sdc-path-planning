# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program
  

## Goal
In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also, the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.


## Solution Summary

To achieve the goal following steps need to be performed for the car:

1. Get car telemetry data from the simulator
2. From the telemetry determine car state:
   * Lane
   * Coordinates
   * Speed
   * Planned path points from the previous planned path
   
3. From telemetry get sensor fusion data - to form predictions - so we can incorporate safety cost to the state machine.

4. Run state machine which: 
   * Determine possible next states
   * For each next state generate possible trajectory
   * Select trajectory with smallest cost.

5. After trajectory is determined send it to the simulator for execution

## Trajectory generation details

Although in the lecture Jerk Minimizing Trajectory approach was used - it was hard to adopt it to this project. (main1.cpp contains the best attempt of now to generate trajectory using this approach). The main issue is that trajectory generated for three parameters coordinates, speed, accelerations - which produce funny behavior for the car. Also having issues to connect previous points to the new model with this approach.   The idea how to solve this it to generate a big set of possible trajectories and select best from them which suites to the previous points. 
However, finally decided to go with spline approach from the walkthrough video. The main idea of the approach is to generate spline ensuring smooth trajectory for the coordinates visited and adjust speed and acceleration to the values required to drive across this trajectory smoothly without jerk. 
Same way is used to generate all trajectories for possible target state - the only difference is target lane number provided to the trajectory generator.

Here is detailed description:
1. Get previous planned points from telemetry and use last two of them - this will ensure that resulting path is smooth and connected.

2. Generate 3 points along desired road lane far away - to ensure smooth lane change or lane driving.  

```cpp
	double target_d = 2.0 + 4.0 * car_lane;
	double s_scale = 50;

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

```

3. Convert points to vehicle coordinates - this just makes math easier as car is driving along X axis.

```cpp
	ConvertToVehicle(pos_x, pos_y, points_x, points_y, angle);
```

4. Using spline library create spline which passes through points from the previous path and currently generated points on the desired lane.
	
```cpp
	spline s;
	s.set_points(points_x, points_y);
```

5. Next we need to identify points along the spline to make sure we are trying to reach target lane speed and preventing huge jerk

```cpp
	// build spline along the points
	spline s;
	s.set_points(points_x, points_y);

	// determine assumed maneuver end coordinates (explicitly actual for lane change)
	double target_x = s_scale;
	double target_y = s(target_x);
	
	// deetermine distance to travel - it is aproximate distance
	double dist = sqrt(target_x * target_x + target_y * target_y);

	double x_start = 0;

	vector<double> points_x_v;
	vector<double> points_y_v;
	
	// for each next point we need to setup speed which is closer to the target speed but do not cause jerk
	double next_speed = target_speed;

	// check if we need to brake to reach target speed
	bool isBrake = (target_speed - speed_end_path) < 0;

	// need to generate 50 points, each next point 0.02 sec away from previous
	for (int i = 0; i < 50; i++)
	{
		// check if acceleration exceed maximum alowed - if so adjust target speed to the value
		// which is possible to reach with maximum allowed acceleration
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

```

6. Convert back to world coordinates so points are ready to be analyzed by cost functions.


```cpp
	spline s;
	s.set_points(points_x, points_y);
```

Please note that such trajectories generated for all possible state transitions. Although we generate 50 points and all used for trajectory cost analysis - not all are sent to the simulator.  We just use couple of new points so together with previous points we have 50 points total.

## Next steps

Although car drives smoothly along highway there are couple of issue enchantments to work with to achive better behaviour.

* Sometimes car stuck in a right lane becouse central is slower and it is costly to change lane there, however left lane is free -so double lane change might be benefitial.   The idea is to make fastest lane as a goal Lane and design additional const fuction to penalize lanes far from goal lane.

* Sometime there is rare issue which make car to change lane to slower lane with car just before it which is extremely hard to debug due to rarity of conditions. Although, this lane changes rarely leads to accident as car often able to adjust speed to avoid accident- the driving behavior looks strange and dangereus. So need to find efficient way to deal with this issue.





