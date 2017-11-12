#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;


// Constant definition
#define LANE_WIDTH (4.0)
#define SIM_STEP (0.02)
#define MIN_SAFE_DIST_FRONT (20)
#define MIN_SAFE_DIST_REAR (20)
#define MAX_DELTA_SPEED (0.45)
#define TARGET_SPEED (49.5)


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

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
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

// we start in lane 1
int ego_lane = 1;
double ref_ego_speed = 0;
double max_s = 6945.554;



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

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          	double ego_x = j[1]["x"];
          	double ego_y = j[1]["y"];
          	double ego_s = j[1]["s"];
          	double ego_d = j[1]["d"];
          	double ego_yaw = j[1]["yaw"];
          	double ego_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];



          	// define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            int prev_path_size = previous_path_x.size();

            //check if we are in the first step
            if (prev_path_size > 0 ){
              ego_s = end_path_s;
            }
            float decel_factor = 0.5;
            bool min_distance_exceeded = false;


            // process sensor fusion information from the other road users
            for (int i = 0; i < sensor_fusion.size(); i++) {

              // ru stands for road user
              int ru_id = sensor_fusion[i][0];
              float ru_vx = sensor_fusion[i][3];
              float ru_vy = sensor_fusion[i][4];
              float ru_speed = sqrt(ru_vx*ru_vx + ru_vy *ru_vy);
              float ru_s = sensor_fusion[i][5];
              float ru_d = sensor_fusion[i][6];

              // check in which lane is the other road user
              int ru_lane = (int) floor(ru_d/LANE_WIDTH);

              if (ru_lane == ego_lane) {

                // estimate car position
                ru_s = ru_s + prev_path_size * SIM_STEP * ru_speed;

                // car is in our lane and closer than the min safe distance, therefore we will slow down
                if ((ru_s - ego_s)<MIN_SAFE_DIST_FRONT && ru_s > ego_s) {

                  min_distance_exceeded = true;

                  // simple P-ish controller for deceleration
                  float delta_speed = ego_speed - ru_speed * 2.24;
                  decel_factor = delta_speed/25 + 0.05; // if delta speed is >50mph, apply full braking

                  //std::cout << "delta speed=" << delta_speed << " decel_factor=" << decel_factor << " ego_speed="<< ego_speed << " target_speed=" << ru_speed * 2.24 << '\n';

                  if ((ru_s - ego_s)< MIN_SAFE_DIST_FRONT/2) {
                    decel_factor = decel_factor * 1.5;
                  }

                  if (decel_factor > 1) {
                    decel_factor = 1;
                  }
                  else if (decel_factor < 0) {
                    decel_factor = 0;
                  }
                }
              }


            }

            /// Lane change
            // only try to do it if there is a vehicle in front and we are not too close to the end of the simulator track (avoid sim errors)
            if (min_distance_exceeded && (ego_s < (max_s - (MIN_SAFE_DIST_FRONT * 2))))
            {
                vector<int> available_lanes;
                // we can only change lanes to adjacent lanes
                if (ego_lane > 0)
                    available_lanes.push_back(ego_lane-1);
                if (ego_lane < 2)
                    available_lanes.push_back(ego_lane+1);

                bool perform_lane_change=false;

                // for every available_lanes and if we are not changing lane
                for ( int j = 0 ; j < available_lanes.size() && !perform_lane_change; j++) {
                    int target_lane = available_lanes[j];

                    double min_distance_forward= max_s;
                    double min_distance_behind= max_s;
                    int closest_ru_behind_i = -1;
                    int closest_ru_front_i = -1;
                    int closest_ru_behind_speed = -1;
                    int closest_ru_front_speed = -1;

                    // process all the sensor fusion objects and find the closest objects in the front and the rear in the target lane
                    for (int i = 0 ; i < sensor_fusion.size(); i++)
                    {
                        int ru_id = sensor_fusion[i][0];
                        float ru_vx = sensor_fusion[i][3];
                        float ru_vy = sensor_fusion[i][4];
                        float ru_speed = sqrt(ru_vx*ru_vx + ru_vy *ru_vy);
                        float ru_s = sensor_fusion[i][5];
                        float ru_d = sensor_fusion[i][6];

                        // check in which lane is the other road user
                        int ru_lane = (int) floor(ru_d/LANE_WIDTH);

                        // estimate the other road user position
                        ru_s = ru_s + prev_path_size * SIM_STEP * ru_speed;

                        if (ru_lane == target_lane)
                        {
                            double delta_dist = ego_s - ru_s;
                            if (delta_dist > 0 ){ // car is behind.
                                if (delta_dist  < min_distance_behind ){
                                    min_distance_behind = delta_dist;
                                    closest_ru_behind_i = i;
                                    closest_ru_behind_speed = ru_speed * 2.24;
                                  }
                            }

                            if (delta_dist < 0 ){
                                if ((delta_dist * -1)< min_distance_forward){
                                    min_distance_forward = delta_dist * -1;
                                    closest_ru_front_i= i;
                                    closest_ru_front_speed = ru_speed * 2.24;
                                  }
                            }
                            if (delta_dist == 0 ){
                                min_distance_behind = min_distance_forward = 0;
                            }
                        }
                    }


                    // check if the lane change conditions are met
                    if (min_distance_behind > MIN_SAFE_DIST_REAR && min_distance_forward > MIN_SAFE_DIST_FRONT){
                        ego_lane = target_lane;
                        perform_lane_change = true;
                    }
                    else if (min_distance_behind > MIN_SAFE_DIST_REAR/3 && min_distance_forward > MIN_SAFE_DIST_FRONT && closest_ru_behind_speed > 0 && closest_ru_behind_speed < (ego_speed + 5) ){
                        ego_lane = target_lane;
                        perform_lane_change = true;
                    }
                    else {
                        cout << "Lane change not possible, too much traffic" << endl;
                    }
                }
            }


            if (min_distance_exceeded) {
              ref_ego_speed = ref_ego_speed - MAX_DELTA_SPEED * decel_factor;
            }
            else if (ref_ego_speed < TARGET_SPEED) {

              // Simple P-ish controller for speed
              float delta_speed = TARGET_SPEED - ego_speed;
              float speed_P_factor = delta_speed/TARGET_SPEED;

              ref_ego_speed = ref_ego_speed + MAX_DELTA_SPEED * speed_P_factor;

              //make sure that we do not exceed the speed limits
              if (ref_ego_speed > TARGET_SPEED) {
                ref_ego_speed = TARGET_SPEED;
              }
            }

            //Create the path points
            double ref_x = ego_x;
            double ref_y = ego_y;
            double ref_yaw = ego_yaw;

            vector<double> pointsX;
            vector<double> pointsY;

            // initialization
            // we consider previous points to smooth the path
            if (prev_path_size <2) {
              double prev_ego_x = ego_x - cos(ego_yaw);
              double prev_ego_y = ego_y - sin(ego_yaw);

              pointsX.push_back(prev_ego_x);
              pointsX.push_back(ego_x);

              pointsY.push_back(prev_ego_y);
              pointsY.push_back(ego_y);
            }
            else {
              ref_x = previous_path_x[prev_path_size-1];
              ref_y = previous_path_y[prev_path_size-1];

              double ref_x_prev = previous_path_x[prev_path_size-2];
              double ref_y_prev = previous_path_y[prev_path_size-2];
              ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

              // to ensure smoothness, we will use these points to get a tangent path to the previous one
              pointsX.push_back(ref_x_prev);
              pointsX.push_back(ref_x);

              pointsY.push_back(ref_y_prev);
              pointsY.push_back(ref_y);
            }

            // using Frenet coordinates we add equally spaced points after the the reference points from the previous trajectory
            // get d position using the ego lane info. Each lane is 4m wide and we add an extra 2m to get the center the lane
            double d_ego_pos = 2.0 + (4.0 * (double)(ego_lane));

            // create points far away that we will use to fit a spline and calculate the actual path
            vector<double> next_wp0 = getXY(ego_s+1*MIN_SAFE_DIST_FRONT*1.5, d_ego_pos, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(ego_s+2*MIN_SAFE_DIST_FRONT*1.5, d_ego_pos, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(ego_s+3*MIN_SAFE_DIST_FRONT*1.5, d_ego_pos, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            pointsX.push_back(next_wp0[0]);
            pointsX.push_back(next_wp1[0]);
            pointsX.push_back(next_wp2[0]);

            pointsY.push_back(next_wp0[1]);
            pointsY.push_back(next_wp1[1]);
            pointsY.push_back(next_wp2[1]);

            // we change the car reference system to a local one with the car looking at an angle of 0deg
            for ( int i = 0 ; i < pointsX.size(); i++)
            {
              // translation + rotation transformation
              double shift_x = pointsX[i] - ref_x;
              double shift_y = pointsY[i] - ref_y;

              pointsX[i] = (shift_x * cos(0 - ref_yaw) - shift_y*sin(0 - ref_yaw));
              pointsY[i] = (shift_x * sin(0 - ref_yaw) + shift_y*cos(0 - ref_yaw));
            }

            // using the spline.h library we create a spline between the reference points and waypoints stored in pointsX and pointsY
            tk::spline s;
            s.set_points(pointsX,pointsY);

            // store the waypoints
            vector<double> next_x_vals;
            vector<double> next_y_vals;

            // we take the previous calculated points that have not been used by the simulator as a starting state for the next path
            for ( int i = 0 ; i < previous_path_x.size(); i++)
            {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            // since speed is dictated by the distance between points in the path, we need to calculate said distance based on the reference speed
            double target_x = 25.0;
            double target_y = s(target_x);
            double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

            double x_increment = 0;

            // we fill the list of path points back again
            for (int i =0 ; i <= 50 - previous_path_x.size();i++)
            {
              double N = (target_dist / (SIM_STEP * ref_ego_speed / 2.24)); // transform mph to m/s -> 2.24 mph to m/s
              double x_point = x_increment+(target_x)/N;
              double y_point = s(x_point);

              x_increment = x_point;

              double x_ref  = x_point;
              double y_ref = y_point;

              // undo the local coordinates transformation back to global coordinates
              // rotation + translation
              x_point = ( x_ref * cos(ref_yaw)-y_ref *sin(ref_yaw));
              y_point = ( x_ref * sin(ref_yaw)+y_ref *cos(ref_yaw));

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);
            }


            json msgJson;

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    cout << "Connected!!!" << endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    cout << "Disconnected" << endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    cout << "Listening to port " << port << endl;
  } else {
    cerr << "Failed to listen to port" << endl;
    return -1;
  }
  h.run();
}
