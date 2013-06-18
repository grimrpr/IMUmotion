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
	static const float time_delta = 0.03;

	measurement.setOnes();
	state.setOnes();

	StateModel<float, state_dimension> state_model;
	state_model.integrate(state, time_delta);
	
	// TEST MeasurementModel

	MeasurementModel<float, state_dimension, measurement_dimension> measurement_model;
	//measurement_model.transformToMeasurement(state, &measurement);

	// TEST UKF
	// testing of the current filter implementation

	UKF<float, state_dimension, measurement_dimension> filter;

	filter.predict(time_delta);
	filter.update(measurement);

	return 0;
}

