/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/


#ifndef MM_H
#define MM_H

#include "Eigen/Core"

template<typename NumType,
	unsigned int state_dimension,
	unsigned int measurement_dimension,
	unsigned int number_of_models=1>
class MeasurementModel
{
public:
	//TODO add sensible confidence value initializiation
	MeasurementModel(){}

	//TODO add real measurement model transformation
	// This function transforms each state_vector into a matching
	// measurement_vector.
	// Implement real measurement model transformation here.
	void transformToMeasurement(const Eigen::Matrix< NumType,
			state_dimension, 2*state_dimension> & state_matrix,
			Eigen::Matrix< NumType, 
			measurement_dimension, 2*state_dimension > *measurement_matrix,
			const unsigned int model_index = 0)
	{
		switch(model_index)
		{
			case 0:
				break;
				//TODO add real measurement model transformation
			default:
				break;
				//measurement_vector = state_vector.tail(measurement_dimension);

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
	
private:

	//this parameter reflects the current confidence in this model and will
	//affect the weighting of the model in the overall orientation calculation
	//the value is inverted and multiplied with the measurement_covariance
	NumType confidence_[number_of_models];
	/* data */
};

/* MM_H */
#endif

