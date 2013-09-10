/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/


#include <iostream>
#include <string>
#include <sstream>
//#include <time>

#include <boost/asio.hpp>

#include "Eigen/Core"
#include "Eigen/Cholesky"

#include "ukf.h"
#include "measurement_model.h"
#include "state_model.h"



int main(int argc, char *argv[])
{

	// TEST of filter
	static const unsigned int state_dimension = 7;
	static const unsigned int measurement_dimension = 9;
	Eigen::Matrix<float, measurement_dimension, 1> measurement;
	Eigen::Matrix<float, state_dimension, 1> state;
	static const float time_delta = 0.01;



	//measurement.setOnes();
	//measurement.x() = 0;
	//measurement.y() = 0;
	//measurement.z() = 1;
	
	measurement(0) = 0;
	measurement(1) = 0;
	measurement(2) = 1;
	measurement(3) = 1;
	measurement(4) = 0;
	measurement(5) = 0;
	measurement(6) = 0;
	measurement(7) = 0;
	measurement(8) = 0;

	state.setZero();
	state(0) = 1;
	state.tail(state_dimension - 4).setZero();

	StateModel<float, state_dimension> state_model;
	std::cout << "state initial: " << std::endl << state << std::endl;
	state_model.integrate(state, time_delta, &state);
	std::cout << "state integrated: " << std::endl << state << std::endl;
	
	// TEST MeasurementModel

	MeasurementModel<float, state_dimension, measurement_dimension> measurement_model;
	//measurement_model.transformToMeasurement(state, &measurement);

	// TEST UKF
	// testing of the current filter implementation

	UKF<float, state_dimension, measurement_dimension> filter;

	Eigen::Matrix<float, state_dimension, 1> initial_state;
	initial_state << 1,0,0,0,1,1,1;

	Eigen::Matrix<float, state_dimension - 1, state_dimension - 1> initial_state_covariance;
	initial_state_covariance.setIdentity();

	Eigen::Matrix<float, state_dimension - 1, state_dimension - 1> process_noise;
	process_noise.setIdentity();

	Eigen::Matrix<float, measurement_dimension, measurement_dimension> measurement_noise;
	measurement_noise.setIdentity();

	filter.Initialize(
			initial_state,
			initial_state_covariance,
			process_noise,
			&measurement_noise);

	//setup serial port
	boost::asio::io_service io_service;
	boost::asio::serial_port port(io_service, "/dev/ttyUSB0");
	port.set_option( boost::asio::serial_port_base::baud_rate(115200) );
	boost::asio::streambuf serial_buffer;

	bool magnet_updated, acceleration_updated, gyro_updated;
	magnet_updated = acceleration_updated = gyro_updated = false;

	//TODO get time delta

	Eigen::Matrix< float, 3, 1> parsed_values;
	std::string line;
	std::string dummy_str, dummy_str2;
	for (;;)
	{
		boost::asio::read_until(port, serial_buffer, '\n');
		std::istream in_stream(&serial_buffer);
		std::getline(in_stream, line);
		std::stringstream ss(line);

		if(line.find("GYRO") != std::string::npos)
		{
			// TODO is this in °/s ?
			// parse values into last three entries of measurement vector
			ss >> dummy_str >> dummy_str2 
				>> parsed_values.x()
			       	>> parsed_values.y()
				>> parsed_values.z();
			measurement.tail<3>() = parsed_values;
			gyro_updated = true;
		}

		else if(line.find("MAG") != std::string::npos)
		{
			// parse values into mid three values of measurement vector
			ss >> dummy_str >> dummy_str2 
				>> parsed_values.x()
			       	>> parsed_values.y()
				>> parsed_values.z();
			measurement.segment<3>(3) = parsed_values.normalized();
			magnet_updated = true;
		}

		else if(line.find("ACC") != std::string::npos)
		{
			// parse values into first three values of measurement vector
			ss >> dummy_str >> dummy_str2 
				>> parsed_values.x()
			       	>> parsed_values.y()
				>> parsed_values.z();
			measurement.head<3>() = parsed_values.normalized();
			acceleration_updated = true;
		}

		if(gyro_updated && acceleration_updated && magnet_updated)
		{
			filter.predict(time_delta);
			filter.update(measurement);

			std::cout << "measurement: " << std::endl << measurement << std::endl;
			magnet_updated = acceleration_updated = gyro_updated = false;
		}

	}

	port.close();

	return 0;
}

