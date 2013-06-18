
IMUmotion
=========

Code for IMU enabled motion capturing and processing. 

Topics:
	- sensor calibration
	- orientation estimation

The sensor calibration is currently only provided for 3-axis magnetometer sensors. It
is based on ellipsoid parameter estimation.
IMU orientation estimation is performed by an unscented kalman filter (UKF). The goal
is to achieve a reusable filter implementation that fully supports quaternion based
orientation tracking and should be able to run in an embedded environment. 
The implementation makes use of the Eigen linear algebra library. Feel encouraged to copy
improve and comment the code. This still is very much work in progress.

author: Benjamin Aschenbrenner

