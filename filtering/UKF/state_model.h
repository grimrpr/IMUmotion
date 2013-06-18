/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/


#ifndef SM_H
#define SM_H

#include "Eigen/Core"
//#include "Eigen/Quaternion"

template<typename NumType,
	unsigned int state_dimension>
class StateModel
{
public:
	StateModel() {}

	// This function predicts for a given state_vector the state after the time 
	// time_delta has passed. The result is stored into state_vector.
	// TODO implement state transition model here.
	// state model:
	// 	- q -unit quaternion (current orientation)
	//	- w -angular rate of each axis
	void integrate(const Eigen::Matrix< NumType,
			state_dimension, 1> & state_vector, NumType time_delta)
	{
		
	}
	
private:
	/* data */
};

/* SM_H */
#endif 

