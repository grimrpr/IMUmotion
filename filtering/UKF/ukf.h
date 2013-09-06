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
	 //standard UKF only
	 //unsigned int number_of_sigma_points = 2*state_dimension,
	 unsigned int number_of_sigma_points = 2*(state_dimension-1),
	 unsigned int number_of_measurement_models = 1 >
class UKF
{
public:
	UKF ()
	{}

	// This function sets the inital filter parameters
	void Initialize(const Eigen::Matrix<NumType, state_dimension, 1> &initial_state, 
			const Eigen::Matrix<NumType, number_of_sigma_points/2, number_of_sigma_points/2> &initial_state_covariance,
			const Eigen::Matrix<NumType, number_of_sigma_points/2, number_of_sigma_points/2> &process_noise,
			const Eigen::Matrix<NumType, measurement_dimension, measurement_dimension> measurement_noise[])
	{
		posteriori_state_vector_ = initial_state;
		posteriori_state_covariance_[0] = initial_state_covariance;
		process_noise_ = process_noise;
		for (unsigned int i = 0; i < number_of_measurement_models; ++i)
		{
			measurement_noise_[i]  = measurement_noise[i];
		}
	}

	void predict(const NumType time_delta)
	{

		// we use Cholesky Decomposition to find the square root of the
		// posteriori state covariance matrix
		Eigen::Matrix< NumType, 
			state_dimension - 1,
		       	state_dimension - 1>  covariance_square_root = (posteriori_state_covariance_[0] + process_noise_).llt().matrixL();

		const NumType pertubation_factor = sqrt(state_dimension);
		covariance_square_root *= pertubation_factor;

		//DEBUG
		std::cout << "covariance_square_root:" << std::endl << covariance_square_root << std::endl;

		Eigen::Matrix< NumType,
			state_dimension,
			number_of_sigma_points> sigma_points;


		// a seccond set of sigma points based on the a priori state estimate,
		// needed for measurement prediction
		Eigen::Matrix< NumType,
			state_dimension,
			number_of_sigma_points> priori_state_sigma_points;

		Eigen::Matrix< NumType,
			measurement_dimension,
			number_of_sigma_points> predicted_measurements[number_of_measurement_models];

		Eigen::Matrix< NumType,
			measurement_dimension,
			number_of_sigma_points> predicted_measurements_centered[number_of_measurement_models];


		//DEBUG
		std::cout << "posteriori_state_vector_:" << std::endl << posteriori_state_vector_ << std::endl;

		// create new sigma point set 
		// TODO put this in external function
		const Eigen::Quaternion<NumType> posteriori_state_vector_quaternion(
					posteriori_state_vector_(0),
					posteriori_state_vector_(1),
					posteriori_state_vector_(2),
					posteriori_state_vector_(3));
		//DEBUG
		std::cout << " posteriori_state_vector_quaternion:" << std::endl << posteriori_state_vector_quaternion.vec() << std::endl;

		for (unsigned int i = 0 ; i < number_of_sigma_points/2; ++i)
		{
			if(! covariance_square_root.col(i).head(3).isZero() )
			{
				const NumType rotation_angle = covariance_square_root.col(i).head(3).norm();
				const Eigen::Matrix< NumType, 3, 1> rotation_axis = 
					covariance_square_root.col(i).head(3).normalized() * sin(rotation_angle/2.0 *M_PI/180.0f);
				//DEBUG
				//std::cout << "rotation angle:" << std::endl << rotation_angle << std::endl;
				//DEBUG
				//std::cout << "rotation axis:" << std::endl << rotation_axis << std::endl;
				Eigen::Quaternion<NumType> pertubation_quaternion( cos(rotation_angle/2.0 *M_PI/180.0f),
						rotation_axis.x(),
						rotation_axis.y(),
						rotation_axis.z() );

				//DEBUG
				//std::cout << "pertubation_quaternion:" << std::endl << pertubation_quaternion.vec() << std::endl;

				pertubation_quaternion = posteriori_state_vector_quaternion * pertubation_quaternion;

				sigma_points.col(i)(0) = pertubation_quaternion.w();
				sigma_points.col(i)(1) = pertubation_quaternion.x();
				sigma_points.col(i)(2) = pertubation_quaternion.y();
				sigma_points.col(i)(3) = pertubation_quaternion.z();
			}
			else
			{
				sigma_points.col(i).head(4) = posteriori_state_vector_.head(4);
			}

			sigma_points.col(i).tail(3) = posteriori_state_vector_.tail(3) + 
				covariance_square_root.col(i).tail(3);

			Eigen::Matrix< NumType, state_dimension, 1> sigma_point_transformed;
			state_model_.integrate(sigma_points.col(i), time_delta, &sigma_point_transformed);
			sigma_points.col(i) = sigma_point_transformed;

			//DEBUG
			//std::cout << "sigma points: " <<  i << std::endl << sigma_points << std::endl;
		}
		
		for (unsigned int i = number_of_sigma_points/2; i < number_of_sigma_points; ++i)
		{
			if(! covariance_square_root.col(i - number_of_sigma_points/2).head(3).isZero() )
			{
				const NumType rotation_angle = covariance_square_root.col(i - number_of_sigma_points/2).head(3).norm();
				const Eigen::Matrix< NumType, 3, 1> rotation_axis = 
					covariance_square_root.col(i - number_of_sigma_points/2).head(3).normalized() * sin(rotation_angle/2.0 *M_PI/180.0f);
				Eigen::Quaternion<NumType> pertubation_quaternion( cos(rotation_angle/2.0 *M_PI/180.0f),
						rotation_axis.x(),
						rotation_axis.y(),
						rotation_axis.z() );

				//DEBUG
				//std::cout << "pertubation_quaternion:" << std::endl << pertubation_quaternion.vec() << std::endl;

				pertubation_quaternion = posteriori_state_vector_quaternion * pertubation_quaternion;

				sigma_points.col(i)(0) = pertubation_quaternion.w();
				sigma_points.col(i)(1) = pertubation_quaternion.x();
				sigma_points.col(i)(2) = pertubation_quaternion.y();
				sigma_points.col(i)(3) = pertubation_quaternion.z();
			}
			else
			{
				sigma_points.col(i).head(4) = posteriori_state_vector_.head(4);
			}

			sigma_points.col(i).tail(3) = posteriori_state_vector_.tail(3) - covariance_square_root.col(i - number_of_sigma_points/2).tail(3);

			Eigen::Matrix< NumType, state_dimension, 1> sigma_point_transformed;
			state_model_.integrate(sigma_points.col(i), time_delta, &sigma_point_transformed);
			sigma_points.col(i) = sigma_point_transformed;

			//DEBUG
			//std::cout << "sigma points: " <<  i << std::endl << sigma_points << std::endl;
	
		}

		//DEBUG
		std::cout << "sigma point set:" << std::endl << sigma_points << std::endl;

		// mean of sigma points is expected state
		// standard UKf
		// priori_state_estimate_ = sigma_points.rowwise().mean();

		// get mean values of measurement parts in sigma_points
		priori_state_estimate_.tail(measurement_dimension) = sigma_points.rowwise().mean().tail(measurement_dimension);

		// get mean values of quaternion parts in sigma_points 
		// TODO find better method for averaging quaternions?
		// TODO put this in external function
		const Eigen::Matrix<NumType, 4, 1> first_pertubated_state_quaternion =
			sigma_points.col(0).head(4);

		Eigen::Matrix<NumType, 4, 1> sum_of_quaternions(0,0,0,0);
		for(unsigned int i = 0; i < number_of_sigma_points; ++i)
		{
			Eigen::Matrix<NumType, 4, 1> pertubated_state_quaternion = 
				sigma_points.col(i).head(4);

			// since a rotation of quaternion q is identical to a rotatiotion
			// of quaternion -q we have to flip the current rotation if
			// it is opposed to the average summation direction
			if(first_pertubated_state_quaternion.dot(pertubated_state_quaternion)
					< 0.0)
			{
				pertubated_state_quaternion = -pertubated_state_quaternion;
			}
			sum_of_quaternions += pertubated_state_quaternion;
		}
		sum_of_quaternions /= static_cast<NumType>(number_of_sigma_points);
		//scale back to unit quaternion
		sum_of_quaternions.normalize();

		priori_state_estimate_.head(4) = sum_of_quaternions;

		//DEBUG
		std::cout << "a priori state vector:" << std::endl << priori_state_estimate_ << std::endl;

		const NumType weight = 1.0 / static_cast<NumType>(number_of_sigma_points);

		// compute covariance of predicted a priori state
		// in the case of normalized quaternions the covariance matrix has
		// dimensions 6x6
		Eigen::Quaternion<NumType> priori_state_quaternion(
				sum_of_quaternions(0),
				sum_of_quaternions(1),
				sum_of_quaternions(2),
				sum_of_quaternions(3));

		Eigen::Matrix<NumType, state_dimension - 1, number_of_sigma_points> independent_sigma_points;

		for (unsigned int i = 0; i < number_of_sigma_points; ++i)
		{
			const Eigen::Quaternion<NumType> sigma_point_quaternion( 
					sigma_points.col(i)(0),
					sigma_points.col(i)(1),
					sigma_points.col(i)(2),
					sigma_points.col(i)(3));
			const Eigen::Quaternion<NumType> error_quaternion = 
				sigma_point_quaternion * priori_state_quaternion.conjugate();

			//DEBUG
			//std::cout << "error_quaternion norm:" << std::endl << error_quaternion.norm() << std::endl;

			independent_sigma_points.col(i).head(3) = ( 2.0*(acos(error_quaternion.w()) * 180.0f/M_PI) ) * error_quaternion.vec();
			independent_sigma_points.col(i).tail(measurement_dimension) = 
				sigma_points.col(i).tail(measurement_dimension) - priori_state_estimate_.tail(measurement_dimension);
		}

		//DEBUG
		std::cout << "independent_sigma_points:" << std::endl << independent_sigma_points << std::endl;

		// only standard UKF
		//Eigen::Matrix< NumType,
		//	state_dimension,
		//	number_of_sigma_points> sigma_points_centered;
		// sigma_points_centered = sigma_points.colwise() - priori_state_estimate_;
		// priori_state_covariance_ = sigma_points_centered
		// 	* sigma_points_centered.transpose()
		// 	* weight
		// 	+ process_noise_;

		priori_state_covariance_ = independent_sigma_points * independent_sigma_points.transpose() * weight;
		//DEBUG
		std::cout << "a priori state vector covariance:" << std::endl << priori_state_covariance_ << std::endl;

		// find sigma points of priori state with priori state covariance
		//covariance_square_root = priori_state_covariance_.llt().matrixL();
		//DEBUG
		//std::cout << "a priori state vector covariance square root:" << std::endl << covariance_square_root << std::endl;

		// TODO do we have to find new pertubated points around estimated state ?
		// //find new set of sigma points from a priori state and priori state covariance
		//for (unsigned int i = 0 ; i < number_of_sigma_points/2; ++i)
		//{
		//	priori_state_sigma_points.col(i) = priori_state_estimate_ +
		//	       (pertubation_factor * covariance_square_root.col(i));
		//		
		//}

		//for (unsigned int i = number_of_sigma_points/2; i < number_of_sigma_points; ++i)
		//{
		//	priori_state_sigma_points.col(i) = priori_state_estimate_ - 
		//		(pertubation_factor * 
		//		 covariance_square_root.col(i - number_of_sigma_points/2));

		//}
		//
		
		priori_state_sigma_points = sigma_points;

		// calculate for each measurement model the following transformations
		for (unsigned int model_index = 0; model_index < number_of_measurement_models; ++model_index)
		{
			// from the state space into the corresponding measurement vector
			measurement_model_.transformToMeasurement(priori_state_sigma_points, 
					model_index, &predicted_measurements[model_index]);

			//DEBUG
			std::cout << "predicted_measurements:" << std::endl << predicted_measurements[model_index] << std::endl;

			// perform state and covariance update for each measurement model
			// mean value is expected measurement
			measurement_estimate_[model_index]  = 
				predicted_measurements[model_index].rowwise().mean();

			//DEBUG
			std::cout << "measurement_estimate_:" << std::endl << measurement_estimate_[model_index] << std::endl;

			// covariance of predicted measurement
			// first get current confidence in model
			const NumType model_confidence = measurement_model_.getModelConfidence(model_index);

			//DEBUG
			std::cout << "model_confidence:" << std::endl << model_confidence << std::endl;

			predicted_measurements_centered[model_index] = 
				predicted_measurements[model_index].colwise() 
				- measurement_estimate_[model_index];

			measurement_covariance_[model_index] = predicted_measurements_centered[model_index]
			       	* predicted_measurements_centered[model_index].transpose()
				* weight + measurement_noise_[model_index] 
				* 1/model_confidence;

			//DEBUG
			std::cout << "measurement_covariance_:" << std::endl << measurement_covariance_[model_index] << std::endl;

			// cross correlation
			cross_covariance_[model_index] = 
				independent_sigma_points * predicted_measurements_centered[model_index].transpose()
				* weight;
			
			//DEBUG
			std::cout << "cross correlation:" << std::endl << cross_covariance_[model_index] << std::endl;

			// Kalman gain
			kalman_gain_[model_index] = cross_covariance_[model_index]
				* (measurement_covariance_[model_index].inverse());
			
			//DEBUG
			std::cout << "kalman_gain_:" << std::endl << kalman_gain_[model_index] << std::endl;

			// update state covariance
			posteriori_state_covariance_[model_index] = 
				priori_state_covariance_ - ( kalman_gain_[model_index] 
						* (cross_covariance_[model_index].transpose()) );

			//DEBUG
			std::cout << "posteriori_state_covariance_:" << std::endl << posteriori_state_covariance_[model_index] << std::endl;
		}

	}

	void update( const Eigen::Matrix< NumType, measurement_dimension, 1> & measurement)
	{
		//a posteriori state estimate
		for ( unsigned int model_index = 0;  model_index < number_of_measurement_models; ++model_index)
		{
			// TODO build in correct weighting of estimated states to gain a single posteriori state
			
			const Eigen::Matrix<NumType, state_dimension - 1, 1> independent_posteriori_state = 
				kalman_gain_[model_index] * (measurement - measurement_estimate_[model_index]);

			//DEBUG
			std::cout << "independent_posteriori_state:" << std::endl << independent_posteriori_state << std::endl;

			const NumType rotation_angle = independent_posteriori_state.head(3).norm();

			//DEBUG
			std::cout << "rotation_angle:" << std::endl << rotation_angle << std::endl;

			const Eigen::Matrix< NumType, 3, 1> rotation_axis = 
				independent_posteriori_state.head(3).normalized() * sin(rotation_angle/2.0 *M_PI/180.0f);

			// add a priori estimate to residual
			posteriori_state_vector_(0) = cos(rotation_angle/2.0);
			posteriori_state_vector_(1) = rotation_axis.x();
			posteriori_state_vector_(2) = rotation_axis.y();
			posteriori_state_vector_(3) = rotation_axis.z();
			posteriori_state_vector_.head(4) += priori_state_estimate_.head(4);
			posteriori_state_vector_.head(4).normalize();
			posteriori_state_vector_.tail(measurement_dimension) = 
				priori_state_estimate_.tail(measurement_dimension) + independent_posteriori_state.tail(measurement_dimension);

			//DEBUG
			std::cout << "posteriori_state_vector_:" << std::endl << posteriori_state_vector_ << std::endl;
		}
	}

private:
//	void transformSigmaPoints();
//	void estimateState();
//	void transformSigmaPoints();
//	void sigmaPointsToMeasurements();
	
	Eigen::Matrix< NumType,
	       	state_dimension,
	       	1> posteriori_state_vector_;

	Eigen::Matrix< NumType,
	       	state_dimension,
	       	1> priori_state_estimate_;

	// only valid in standard UKF
	// Eigen::Matrix< NumType, 
	// 	state_dimension, 
	// 	state_dimension> posteriori_state_covariance_[number_of_measurement_models];
	
	// since our quaternions are normalized in the state vector the four quaternion elements are
	// stripped of one dimension which reflects in the covariance and noise
	// representation of the state
	Eigen::Matrix< NumType, 
		state_dimension - 1, 
		state_dimension - 1> posteriori_state_covariance_[number_of_measurement_models];

	// only valid in standard UKF
	// Eigen::Matrix< NumType, 
	// 	state_dimension, 
	// 	state_dimension> priori_state_covariance_;

	Eigen::Matrix< NumType, 
		state_dimension - 1, 
		state_dimension - 1> priori_state_covariance_;

	// only valid in standard UKF
	// Eigen::Matrix< NumType,
	//        	state_dimension,
	//        	state_dimension> process_noise_;
	

	Eigen::Matrix< NumType,
	       	state_dimension - 1,
	       	state_dimension - 1> process_noise_;

	Eigen::Matrix< NumType,
	       	measurement_dimension,
	       	1> measurement_estimate_[number_of_measurement_models];

	Eigen::Matrix< NumType, 
		measurement_dimension, 
		measurement_dimension> measurement_covariance_[number_of_measurement_models];

	Eigen::Matrix< NumType, 
		measurement_dimension, 
		measurement_dimension> measurement_noise_[number_of_measurement_models];

	// only valid in standard UKF
	// Eigen::Matrix< NumType,
	// 	state_dimension,
	// 	measurement_dimension> cross_covariance_[number_of_measurement_models];

	Eigen::Matrix< NumType,
		state_dimension - 1,
		measurement_dimension> cross_covariance_[number_of_measurement_models];

	// only valid in standard UKF
	// Eigen::Matrix< NumType,
	// 	state_dimension,
	// 	measurement_dimension> kalman_gain_[number_of_measurement_models];

	Eigen::Matrix< NumType,
		state_dimension - 1,
		measurement_dimension> kalman_gain_[number_of_measurement_models];

	MeasurementModel<NumType,
	       	state_dimension, 
		measurement_dimension,
	       	number_of_sigma_points,	
		number_of_measurement_models> measurement_model_;

	StateModel<NumType, state_dimension> state_model_;


};

/* UKF_H */
#endif 
