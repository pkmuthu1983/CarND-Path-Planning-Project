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

#define PRINTSF
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

void printpoints(vector<double> x, vector<double> y, string name) {
#ifdef DEBUG
	cout << "Points " << name << endl;
	for (int i=0; i<x.size(); i++){
		cout << "idx - " << i << ", (" << x[i] << ", " << y[i] << ")" << endl;
	}
#endif
}

void printfusiondata(vector< vector<double> > sf, int iter) {
#ifdef DEBUG
	cout << "Iter " << iter << ", Fusion data size" << sf.size() << endl;
	for (int i=0; i<sf.size(); i++) {
		for (int j=0; j<7; j++) {
			cout << sf[i][j] << ",";
		}
		cout << endl;
	}
#endif
}

void printlanedata(int curlane, vector<int> nextlanes,
    vector<int> lanegd, int bestlane, vector<vector<double> > sf, int iter) {

#ifdef DEBUG
	assert (nextlanes.size() == lanegd.size());

	cout << "Considering lane change. " << (bestlane == -1 ? "None found": "Good lane found") << endl;
	cout << "cur lane " << curlane << endl;

	for (int i=0; i<lanegd.size(); i++) {
		cout << nextlanes[i] << ": " << lanegd[i] << ",";
	}
	cout << "bestlane " << bestlane << endl;

	cout << "Iter " << iter << ", Fusion data size" << sf.size() << endl;

	for (int i=0; i<sf.size(); i++) {
		for (int j=0; j<7; j++) {
			cout << sf[i][j] << ",";
		}
		cout << endl;
	}
#endif
}

int findlane(double d) {
	int lane = 0;
	if (d >= 0 && d<4) {
		lane = 0;	
	} else if (d >= 4 && d < 8) {
		lane = 1;
	} else if (d >= 8 && d < 12) {
		lane = 2;
	} else {
		lane = -1;
	}
	return lane;
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
  double ref_v = 0;
  int curlane = -1, nextlane = -1, lane_change_ongoing = 0, last_lane_change = 0;
  int iter = 0;

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

  h.onMessage([&iter, &ref_v, &curlane, &nextlane, &lane_change_ongoing, &last_lane_change, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;

    iter++;

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

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

		vector<double> pts_x;
		vector<double> pts_y;

		double dist_inc = 0.44;

		double ref_x;
		double ref_y;
		double angle;

		int path_size = previous_path_x.size();

		if (curlane == -1) {
			curlane = findlane(car_d);
		}

		double car_s_orig = car_s;

		if (path_size > 0) {
			car_s = end_path_s;
		}

		bool tooclose = false;
		double car_in_front_s = 1000;
		int car_in_front_id = -1;
		double car_in_front_speed = 0;

		for (int i=0; i< sensor_fusion.size(); i++){
			double d = sensor_fusion[i][6];
			if ((d >= car_d - 2 && d <= car_d + 2) ||
			    (findlane(car_d) == findlane(d)) ) {
				double vx = sensor_fusion[i][3];
				double vy = sensor_fusion[i][4];
				double opp_car_speed = sqrt(vx*vx + vy*vy);
				double opp_car_s = sensor_fusion[i][5];
				if (opp_car_s > car_s && opp_car_s - car_s < 30) {
					tooclose = true;
					if (opp_car_s < car_in_front_s) {
						car_in_front_s = opp_car_s;
						car_in_front_id = i;
					}
				}
			}
		}

		if (car_in_front_id != -1) {
			double vx = sensor_fusion[car_in_front_id][3];
			double vy = sensor_fusion[car_in_front_id][4];
			car_in_front_speed = sqrt(vx*vx + vy*vy) * 2.23 - 2;
		}

		if (tooclose) {
			if (car_in_front_id != -1 && ref_v < car_in_front_speed) {
				ref_v = car_in_front_speed;
#ifdef DEBUG	
				cout << "Following car " << car_in_front_id << " with " << car_in_front_speed << endl;
#endif
			} else {
				ref_v -= 0.5;
			}
		} else if (ref_v < 49) {
			ref_v += 0.5;
		}

		if (tooclose && !lane_change_ongoing && (iter - last_lane_change > 250) ) {
			vector<int> nextlanes;
			if (curlane == 0) {
				nextlanes.push_back(1);
			} else if (curlane == 1) {
				nextlanes.push_back(0);
				nextlanes.push_back(2);
			} else if (curlane == 2) {
				nextlanes.push_back(1);
			} else {
				assert(false);
			}

			vector<int> lanegood;
	
			int bestlane = -1;

			for (int j = 0; j < nextlanes.size(); j++) {
				int thislane = nextlanes[j];
				lanegood.push_back(1);

				for (int i=0; i< sensor_fusion.size(); i++){
					double d = sensor_fusion[i][6];
					if (d >= (2 + 4 * thislane - 2) && d <= (2 + 4 * thislane + 2)) {
						double vx = sensor_fusion[i][3];
						double vy = sensor_fusion[i][4];
						double opp_car_speed = sqrt(vx*vx + vy*vy);
						double opp_car_s = sensor_fusion[i][5];
						opp_car_s += opp_car_speed * 0.02 * path_size;

						if (abs(opp_car_s - car_s) < 30 ||
						    abs(opp_car_s - car_s_orig) < 30) {
							lanegood[j] = 0;
							break;
						}
					}
				}
			}

			for (int j = 0; j < nextlanes.size(); j++) {
				if (lanegood[j] == 1) {
					bestlane = nextlanes[j];
					break;
				}
			}

			if (bestlane != -1) {
				lane_change_ongoing = 1;
				nextlane = bestlane;
				last_lane_change = iter;
			}

			printlanedata(curlane, nextlanes, lanegood, bestlane, sensor_fusion, iter);
		}

		double lane = -1;

		if (lane_change_ongoing == 1) {
			lane = nextlane;
		} else {
			lane = curlane;
		}

		if (nextlane != -1 && curlane != nextlane) {
			cout << "Lane changing from " << curlane << " to " << nextlane << ". car_s " << car_s << ", car_d " << car_d << endl;
			if (findlane(car_d) == nextlane) {
				curlane = nextlane;
				lane_change_ongoing = 0;
				nextlane = -1;
#ifdef DEBUG	
				cout << "Lane change complete" << endl;
#endif
			}
		}


		double ref_speed_ms = ref_v/2.23;

		if(path_size < 2) {
			angle = deg2rad(car_yaw);
			double prv_x = car_x - dist_inc * cos(angle);
			double prv_y = car_y - dist_inc * sin(angle);
			pts_x.push_back(prv_x);
			pts_x.push_back(car_x);
			pts_y.push_back(prv_y);
			pts_y.push_back(car_y);
			ref_x = car_x;
			ref_y = car_y;
		} else {
			double ref_x_1, ref_x_2, ref_y_1, ref_y_2;
			ref_x_2 = previous_path_x[path_size-2];
			ref_x_1 = previous_path_x[path_size-1];
			ref_y_2 = previous_path_y[path_size-2];
			ref_y_1 = previous_path_y[path_size-1];

			pts_x.push_back(ref_x_2);
			pts_x.push_back(ref_x_1);
			pts_y.push_back(ref_y_2);
			pts_y.push_back(ref_y_1);
			angle = atan2(ref_y_1 - ref_y_2, ref_x_1 - ref_x_2);
			ref_x = ref_x_1;
			ref_y = ref_y_1;
		}
	
		vector<double> nxt_wp0 = getXY(car_s + 30, 2 + 4 * lane,
		    map_waypoints_s, map_waypoints_x, map_waypoints_y);
		vector<double> nxt_wp1 = getXY(car_s + 60, 2 + 4 * lane,
		    map_waypoints_s, map_waypoints_x, map_waypoints_y);
		vector<double> nxt_wp2 = getXY(car_s + 90, 2 + 4 * lane,
		    map_waypoints_s, map_waypoints_x, map_waypoints_y);

		pts_x.push_back(nxt_wp0[0]);
		pts_x.push_back(nxt_wp1[0]);
		pts_x.push_back(nxt_wp2[0]);

		pts_y.push_back(nxt_wp0[1]);
		pts_y.push_back(nxt_wp1[1]);
		pts_y.push_back(nxt_wp2[1]);

		for (unsigned int i = 0; i < pts_x.size(); i++) {
			double shiftx = pts_x[i] - ref_x;
			double shifty = pts_y[i] - ref_y;
		
			pts_x[i] = (shiftx * cos(0-angle) - shifty * sin(0-angle));	
			pts_y[i] = (shiftx * sin(0-angle) + shifty * cos(0-angle));
		}

		assert (pts_x.size() == pts_y.size());

		tk::spline s;
		s.set_points(pts_x, pts_y);

          	for(int i = 0; i < path_size; i++)
		{
			next_x_vals.push_back(previous_path_x[i]);
			next_y_vals.push_back(previous_path_y[i]);
		}

		double target_x = 45.0;
		double target_y = s(target_x);
		double target_dist = sqrt((target_x*target_x) + (target_y)*(target_y));

		double cur_x = 0, cur_y = 0, cur_x_T = 0, cur_y_T = 0;
		double N = target_dist/(0.02*ref_speed_ms);
		double inc = target_x/N;

		for (int i = 1; i <= 60 - path_size; i++) {
			cur_x += inc;
			cur_y = s(cur_x);
			cur_x_T = ref_x + (cur_x * cos(angle) - cur_y * sin(angle));
			cur_y_T = ref_y + (cur_x * sin(angle) + cur_y * cos(angle));
			next_x_vals.push_back(cur_x_T);
			next_y_vals.push_back(cur_y_T);
		}

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
