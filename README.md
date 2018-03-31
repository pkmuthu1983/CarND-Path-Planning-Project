# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

[video1]: ./path_planning.mp4 "Project Video"

### Source code:
src/main.cpp implements path planner code. It can be run using build/path_planner

### Model Description:

1) The path planner generates jerk minimizing trajectory by fitting a polynomial (via spline library) to
a set of points that include the current car position and future positions that are 30, 60, and 90m
in front of the car. Then based on reference velocity, evenly spaced points are obtained for the car
to follow in future.

2) If the car senses another car in front of it, then the velocity is reduced, but it will try to 
match the speed of the car in front of it. Refer lines 334 to 374 of main.cpp for the implementation.

3) The ego car will try to switch lanes, whenever it finds another in front of it. It will try to 
find the best lane, based on its current lane and the traffic in nearby lanes. Whenever an opportunity
exists, it switches lanes. Moreover, the car is hard-coded NOT to change lanes too frequently (should
be atleast 5 seconds gap). Refer lines 376 to 449, main.cpp.

4) The trajectories for lane change is also obtained using the spline technique mentioned in 1. 
In this case however, the future points are chosen from the desired lane.

### Results:
The car can successfully drive 4.32 miles without any incidents as shown in the video below:
[link to video](./path_planning.mp4)

