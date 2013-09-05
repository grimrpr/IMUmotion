/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/


#ifndef MM_H
#define MM_H

#include "Eigen/Core"
#include "Eigen/Geometry"

enum ModelIdentifier{
	ACCELEROMETER,
	GYROSCOPE,
	MAGNETOMETER
};

template<typename NumType,
	unsigned int state_dimension,
	unsigned int measurement_dimension,
	//only in standard UKF
	unsigned int number_of_sigma_points = 2*state_dimension,
	unsigned int number_of_models=1>
class MeasurementModel
{
public:
	//TODO add sensible confidence value initializiation
	MeasurementModel() : 
	north(1,0,0),
	east(0,1,0),
	down(0,0,1){}

	// This function transforms each state_vector into a matching
	// measurement_vector.
	// Implement real measurement model transformation here.
	void transformToMeasurement(
			const Eigen::Matrix< NumType, state_dimension, number_of_sigma_points> & state_matrix,
			const unsigned int model_index,
			Eigen::Matrix< NumType, measurement_dimension, number_of_sigma_points > *measurement_matrix)
	{
		switch(model_index)
		{
			case ACCELEROMETER:

				for (unsigned int i = 0; i < number_of_sigma_points; ++i)
				{
					const Eigen::Matrix< NumType, state_dimension, 1> state_vector = state_matrix.col(i);
					const Eigen::Quaternion<NumType> state_quaternion( 
						state_vector(0), 
						state_vector(1), 
						state_vector(2), 
						state_vector(3) );

					measurement_matrix->col(i) = state_quaternion._transformVector(down);
				}

				break;

			case MAGNETOMETER:

				for (unsigned int i = 0; i < number_of_sigma_points; ++i)
				{
					const Eigen::Matrix< NumType, state_dimension, 1> state_vector = state_matrix.col(i);
					const Eigen::Quaternion<NumType> state_quaternion( 
						state_vector(0), 
						state_vector(1), 
						state_vector(2), 
						state_vector(3) );

					measurement_matrix->col(i) = state_quaternion._transformVector(north);
				}

				break;

			case GYROSCOPE:
			default:
				measurement_matrix->bottomLeftCorner(measurement_dimension, number_of_sigma_points) = state_matrix.bottomLeftCorner(
						measurement_dimension, 
						number_of_sigma_points );
				break;

		}
	}

	NumType getModelConfidence(const unsigned int model_index)
	{
		if(model_index < number_of_models)
			return confidence_[model_index];
		else
			return 0;
	}

	void setModelConfidence(const unsigned int model_index, NumType new_confidence)
	{
		if(model_index < number_of_models)
			confidence_[model_index] = new_confidence;
	}
	
public:
	const Eigen::Matrix< NumType, 3, 1> north;
	const Eigen::Matrix< NumType, 3, 1> east;
	const Eigen::Matrix< NumType, 3, 1> down;

private:

	//this parameter reflects the current confidence in this model and will
	//affect the weighting of the model in the overall orientation calculation
	//the value is inverted and multiplied with the measurement_covariance
	NumType confidence_[number_of_models];
	/* data */
};

/* MM_H */
#endif

