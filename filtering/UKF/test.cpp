/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/


#include <iostream>

#include "Eigen/Core"
#include "Eigen/Cholesky"

#include "ukf.h"
#include "measurement_model.h"
#include "state_model.h"


int main(int argc, char *argv[])
{

	// TEST Eigen
	// we do some Eigen function testing here
	Eigen::Matrix< float, 3, 3> matrix;
	matrix << 4 , -1, 2, -1, 6, 0, 2, 0, 5;
	
	// Cholesky Decomposition
	Eigen::Matrix< float, 3, 3> matrix_square_root = matrix.llt().matrixL();

	std::cout << "The matrix 1 is:" << std::endl 
		<< matrix << std::endl;
	std::cout << "The matrix 1 square root is:" << std::endl 
		<< matrix_square_root << std::endl;

	// TEST StateModel
	static const unsigned int state_dimension = 7;
	static const unsigned int measurement_dimension = 3;
	Eigen::Matrix<float, measurement_dimension, 1> measurement;
	Eigen::Matrix<float, state_dimension, 1> state;
	//static const float time_delta = 0.03;
	static const float time_delta = 1.0;

	//measurement.setOnes();
	measurement.x() = 0;
	measurement.y() = 0;
	measurement.z() = 1;

	state.setZero();
	state(0) = 1;
	state.tail(measurement_dimension) = Eigen::Vector3f(1,1,1);

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
	initial_state << 1,0,0,0,3,3,3;
	Eigen::Matrix<float, state_dimension - 1, state_dimension - 1> initial_state_covariance;
	initial_state_covariance.setIdentity();
	Eigen::Matrix<float, state_dimension - 1, state_dimension - 1> process_noise;
	process_noise.setZero();
	Eigen::Matrix<float, measurement_dimension, measurement_dimension> measurement_noise;
	measurement_noise.setZero();

	filter.Initialize(
			initial_state,
			initial_state_covariance,
			process_noise,
			&measurement_noise);

	filter.predict(time_delta);
	filter.update(measurement);

	return 0;
}

