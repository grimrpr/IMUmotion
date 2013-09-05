/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/


#ifndef SM_H
#define SM_H

#include "Eigen/Core"
#include "Eigen/Geometry"

template<typename NumType,
	unsigned int state_dimension>
class StateModel
{
public:
	StateModel() {}

	// This function predicts for a given state_vector the state after the time 
	// time_delta has passed. The result is stored into result_state_vector.
	// state model:
	// 	- q -unit quaternion (current orientation)
	//	- w -angular rate of each axis
	void integrate(const Eigen::Matrix< NumType, state_dimension, 1> & state_vector, 
			const NumType time_delta,
			Eigen::Matrix< NumType, state_dimension, 1> * result_state_vector)
	{
		if(!state_vector.tail(3).isZero())
		{
			const Eigen::Quaternion<NumType> state_quaternion( 
					state_vector(0),
					state_vector(1),
					state_vector(2),
					state_vector(3) );
			
			// Create quaternion that rotates the state_vector quaternion into
			// the new position predicted at the current rotation rate and time
			// time_delta.
			const NumType rotation_angle = state_vector.tail(3).norm() * time_delta;
			const Eigen::Matrix< NumType, 3, 1> rotation_axis = 
				state_vector.tail(3).normalized() * sin(rotation_angle/2);
			const Eigen::Quaternion<NumType> rotation_quaternion( cos(rotation_angle/2),
					rotation_axis.x(),
					rotation_axis.y(),
					rotation_axis.z() );

			// TODO do we need to normalize the result?
			// perform rotation
			const Eigen::Quaternion<NumType> result_quaternion = 
				state_quaternion * rotation_quaternion;

			(*result_state_vector)(0) = result_quaternion.w();
			(*result_state_vector)(1) = result_quaternion.x();
			(*result_state_vector)(2) = result_quaternion.y();
			(*result_state_vector)(3) = result_quaternion.z();
			result_state_vector->tail(3) = state_vector.tail(3);
		}
		else
		{
			*result_state_vector = state_vector;
		}
	
	}
	
private:
	/* data */
};

/* SM_H */
#endif 

