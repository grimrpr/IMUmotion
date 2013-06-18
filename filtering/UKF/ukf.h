/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/


#ifndef	UKF_H
#define	UKF_H

#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/Geometry"

#include "measurement_model.h"
#include "state_model.h"

#include <math.h>
#include <iostream>

template < typename NumType, 
	 unsigned int state_dimension, 
	 unsigned int measurement_dimension,
	 unsigned int number_of_measurement_models = 1 >
class UKF
{
public:
	UKF ()
		//TODO initialize measurement noise
		//TODO initialize process noise
	{}

	void predict(const NumType time_delta)
	{

		//TODO correct the initialization of state_vector and covariance
		posteriori_state_vector_ = 
		Eigen::Matrix< NumType, state_dimension, 1>::Zero();
		posteriori_state_covariance_[0].setIdentity();

	//	posteriori_state_covariance_ << 
	//			       	0, 5, 4, -1, 2, -1, 6, 
	//			       	0, 5, 4, -1, 2, -1, 6, 
	//			       	0, 5, 4, -1, 2, -1, 6, 
	//			       	0, 5, 4, -1, 2, -1, 6, 
	//			       	0, 5, 4, -1, 2, -1, 6, 
	//			       	0, 5, 4, -1, 2, -1, 6, 
	//			       	0, 5, 4, -1, 2, -1, 6; 
	

		// Cholesky Decomposition to find matrix S so that
		// posteriori_state_covariance = S * S^T
		// this is also known as the square root of a matrix
		Eigen::Matrix< NumType, 
			state_dimension,
		       	state_dimension >  covariance_square_root = posteriori_state_covariance_[0].llt().matrixL();

		Eigen::Matrix< NumType,
			state_dimension,
			2* state_dimension> sigma_points;
		Eigen::Matrix< NumType,
			state_dimension,
			2* state_dimension> sigma_points_centered;


		// a seccond set of sigma points based on the a priori state estimate,
		// needed for measurement prediction
		Eigen::Matrix< NumType,
			state_dimension,
			2* state_dimension> priori_state_sigma_points;

		Eigen::Matrix< NumType,
			measurement_dimension,
			2* state_dimension> predicted_measurements[number_of_measurement_models];

		Eigen::Matrix< NumType,
			measurement_dimension,
			2* state_dimension> predicted_measurements_centered[number_of_measurement_models];


		//DEBUG
		std::cout << "current estimated state:" << std::endl << posteriori_state_vector_ << std::endl;

		// create new sigma point set 
		// TODO put the following constants into constructor list
		const unsigned int number_of_sigma_points = sigma_points.cols(); 
		const NumType pertubation_factor = sqrt(state_dimension);
		const NumType weight = 1.0 / (2 * state_dimension);

		//TODO add state model prediction step
		for (unsigned int i = 0 ; i < number_of_sigma_points/2; ++i)
		{
			sigma_points.col(i) = posteriori_state_vector_ +
			       (pertubation_factor * covariance_square_root.col(i));
			//state_model_.integrate(&sigma_points.col(i), time_delta);
		}
		
		for (unsigned int i = number_of_sigma_points/2; i < number_of_sigma_points; ++i)
		{
			sigma_points.col(i) = posteriori_state_vector_ - 
				(pertubation_factor * 
				 covariance_square_root.col(i - number_of_sigma_points/2));
			//state_model_.integrate(&sigma_points.col(i), time_delta);
		}

		//DEBUG
		std::cout << "sigma point set:" << std::endl << sigma_points << std::endl;

		// mean of sigma points is expected state
		//TODO optimize for quaternions just taking the mean is not realy
		// what we want when dealing with quaternions
		priori_state_estimate_ = sigma_points.rowwise().mean();

		//DEBUG
		std::cout << "a priori state vector:" << std::endl << priori_state_estimate_ << std::endl;

		// compute covariance of predicted state
		sigma_points_centered = sigma_points.colwise() - priori_state_estimate_;
		priori_state_covariance_ = sigma_points_centered
			* sigma_points_centered.transpose()
			* weight
			+ process_noise_;
		//DEBUG
		std::cout << "a priori state vector covariance:" << std::endl << priori_state_covariance_ << std::endl;

		// find sigma points of priori state with priori state covariance
		covariance_square_root = priori_state_covariance_.llt().matrixL();

		//find new set of sigma points from a priori state and priori state covariance
		for (unsigned int i = 0 ; i < number_of_sigma_points/2; ++i)
		{
			priori_state_sigma_points.col(i) = priori_state_estimate_ +
			       (pertubation_factor * covariance_square_root.col(i));
				
		}

		for (unsigned int i = number_of_sigma_points/2; i < number_of_sigma_points; ++i)
		{
			priori_state_sigma_points.col(i) = priori_state_estimate_ - 
				(pertubation_factor * 
				 covariance_square_root.col(i - number_of_sigma_points/2));

		}

		// calculate for each measurement model the following transformations
		for (unsigned int model_index = 0; model_index < number_of_measurement_models; ++model_index)
		{
			// from the state space into the corresponding measurement vector
			measurement_model_.transformToMeasurement(priori_state_sigma_points, 
					&predicted_measurements[model_index], model_index);


			// perform state and covariance update for each measurement model
			// mean value is expected measurement
			measurement_estimate_[model_index]  = 
				predicted_measurements[model_index].rowwise().mean();

			//covariance of predicted measurement
			//first get current confidence in model
			const NumType model_confidence = measurement_model_.getModelConfidence(model_index);

			predicted_measurements_centered[model_index] = 
				predicted_measurements[model_index].colwise() 
				- measurement_estimate_[model_index];

			measurement_covariance_[model_index] = predicted_measurements_centered[model_index]
			       	* predicted_measurements_centered[model_index].transpose()
				* weight + measurement_noise_[model_index] 
				* 1/model_confidence;

			// TODO Check: is it correct that we use priori_state_sigma_points
			// instead of sigma_points_centered ?
			// cross covariance
			cross_covariance_[model_index] = 
				sigma_points_centered * predicted_measurements_centered[model_index].transpose()
				* weight;

			//Kalman gain
			kalman_gain_[model_index] = cross_covariance_[model_index]
				* (measurement_covariance_[model_index].inverse());

			//update state covariance
			// TODO get correct posteriori state covariance
			posteriori_state_covariance_[model_index] = 
				priori_state_covariance_ - ( kalman_gain_[model_index] 
						* (cross_covariance_[model_index].transpose()) );
		}

	}

	void update( const Eigen::Matrix< NumType, measurement_dimension, 1> & measurement)
	{
		//a posteriori state estimate
		for ( unsigned int model_index = 0;  model_index < number_of_measurement_models; ++model_index)
		{
			//TODO build in correct weighting of estimated states to gain
			// one state
			posteriori_state_vector_ = priori_state_estimate_ 
				+ kalman_gain_[model_index] * (measurement - measurement_estimate_[model_index]);
		}
	}

private:
//	void transformSigmaPoints();
//	void estimateState();
//	void transformSigmaPoints();
//	void sigmaPointsToMeasurements();
	
	MeasurementModel<NumType,
	       	state_dimension, 
		measurement_dimension, 
		number_of_measurement_models> measurement_model_;

	StateModel<NumType, state_dimension> state_model_;

	Eigen::Matrix< NumType,
	       	state_dimension,
	       	1> posteriori_state_vector_;

	Eigen::Matrix< NumType,
	       	state_dimension,
	       	1> priori_state_estimate_;

	Eigen::Matrix< NumType, 
		state_dimension, 
		state_dimension> posteriori_state_covariance_[number_of_measurement_models];

	Eigen::Matrix< NumType, 
		state_dimension, 
		state_dimension> priori_state_covariance_;

	Eigen::Matrix< NumType,
	       	state_dimension,
	       	state_dimension> process_noise_;

	Eigen::Matrix< NumType,
	       	measurement_dimension,
	       	1> measurement_estimate_[number_of_measurement_models];

	Eigen::Matrix< NumType, 
		measurement_dimension, 
		measurement_dimension> measurement_covariance_[number_of_measurement_models];

	Eigen::Matrix< NumType, 
		measurement_dimension, 
		measurement_dimension> measurement_noise_[number_of_measurement_models];

	Eigen::Matrix< NumType,
		state_dimension,
		measurement_dimension> cross_covariance_[number_of_measurement_models];

	Eigen::Matrix< NumType,
		state_dimension,
		measurement_dimension> kalman_gain_[number_of_measurement_models];

};

/* UKF_H */
#endif 
