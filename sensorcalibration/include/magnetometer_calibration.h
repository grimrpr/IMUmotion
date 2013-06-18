/*
 * author: Benjamin Aschenbrenner
 * year: 2013
 * Feel free to copy and modify just include this header.
*/

#include <iostream>
#include <random>
#include <fstream>

#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Eigenvalues"

//TODO check parameter type:  use unsigned short or is uint16 available?
extern "C" void make_calibration(unsigned short xyz_values[])
{

  //defines how many samples will be used and stored per calibratin cycle
  static const unsigned int number_of_points = 200;

  //how many samples have been already read 
  static unsigned int number_of_read_points = 0;

  static bool initial_calibration = true;
  static Matrix<float, 9, 1> parameters;
  static Matrix<float, 3, 1> center;
  static Matrix<float, 3, 3> eigenvalues;
  static Matrix<float, 3, 3> eigenvectors;
  static Matrix<float, 3, 3> scale_matrix;
  static Array<float, 3, number_of_points> test_points_ellipsoid;


  Matrix<float, 3, 1> in_point;

  in_point.x() = xyz_values[0];
  in_point.y() = xyz_values[1];
  in_point.z() = xyz_values[2];

  if(number_of_read_points < number_of_points)
  {

	  test_points_ellipsoid.col(number_of_read_points) = in_point;
	  ++number_of_read_points;
  }

  else
  {
	//get new calibration parameters
	number_of_read_points = 0;

  	LinearLeastSquares_QuadraticForm< float, number_of_points>( test_points_ellipsoid, &parameters);
  	GetEllipsoidParameters(parameters, &center, &eigenvalues, &eigenvectors);

  	//remove hard iron and bias
  	test_points_ellipsoid.colwise() -= center.array();

  	eigenvalues(0,0) = internal::sqrt(1.0/eigenvalues(0,0));
  	eigenvalues(1,1) = internal::sqrt(1.0/eigenvalues(1,1));
  	eigenvalues(2,2) = internal::sqrt(1.0/eigenvalues(2,2));
 
  	scale_matrix(0,0) = (1.0/eigenvalues(0,0));
  	scale_matrix(1,1) = (1.0/eigenvalues(1,1));
  	scale_matrix(2,2) = (1.0/eigenvalues(2,2));

	initial_calibration = false;
  }

  if(!initial_calibration)
  {
  	//scale to unitvector
  	in_point = eigenvectors * scale_matrix * eigenvectors.transpose() * in_point;

	xyz_values[0] = in_point.x();
	xyz_values[1] = in_point.y();
	xyz_values[2] = in_point.z();

 // 	std::ofstream file3("points_retransformed.txt");
 // 	if(file3.is_open())
 //	 {
 //	   file3 << test_points_ellipsoid.transpose() << std::endl;
 //	   file3.close();
 //	 }
  }

}

using namespace Eigen;
template < typename ScalarType, unsigned int number_of_points >
void LinearLeastSquares_QuadraticForm(const Array< ScalarType, 3, number_of_points> & data_points, Matrix< ScalarType, 9, 1> * result)
{
  
  Matrix< ScalarType, number_of_points, 9> X;

  //x²
  X.col(0) = data_points.row(0)*data_points.row(0);

  //y²
  X.col(1) = data_points.row(1)*data_points.row(1);

  //z²
  X.col(2) = data_points.row(2)*data_points.row(2);
  
  //2xy
  X.col(3) = 2*data_points.row(0)*data_points.row(1);

  //2xz
  X.col(4) = 2*data_points.row(0)*data_points.row(2);

  //2yz
  X.col(5) = 2*data_points.row(1)*data_points.row(2);

  //2*x
  X.col(6) = 2*data_points.row(0);
  
  //2*y
  X.col(7) = 2*data_points.row(1);
  
  //2*z
  X.col(8) = 2*data_points.row(2);

  //result vector
  Matrix< ScalarType, number_of_points, 1> T = 
    Array< ScalarType, number_of_points, 1>::Constant(1);

  //T = data_points.row(2)*data_points.row(2);
  //std::cout << X << std::endl;

  //least squares
  *result = (X.transpose() * X).inverse() * X.transpose()  * T;
}

template<class ScalarType>
void GetEllipsoidParameters(const Matrix< ScalarType, 9, 1> &lsq_parameters, 
    Matrix< ScalarType, 3, 1> * center,
    Matrix< ScalarType, 3, 3 > * eigenvalue_matrix,
    Matrix< ScalarType, 3, 3 > * eigenvectors )
{
  Matrix< ScalarType, 4, 4> algebraic_form;
  algebraic_form << lsq_parameters(0), lsq_parameters(3), lsq_parameters(4), lsq_parameters(6),
                    lsq_parameters(3), lsq_parameters(1), lsq_parameters(5), lsq_parameters(7),
                    lsq_parameters(4), lsq_parameters(5), lsq_parameters(2), lsq_parameters(8),
                    lsq_parameters(6), lsq_parameters(7), lsq_parameters(8), -1;

  Matrix< ScalarType, 3, 1> xyz_vector;
  xyz_vector << lsq_parameters(6), lsq_parameters(7), lsq_parameters(8);

  *center = (algebraic_form.topLeftCorner(3,3) * (-1) ).inverse() * xyz_vector;

  Matrix< ScalarType, 4, 4> translation_matrix = Matrix< ScalarType, 4, 4>::Identity();
  translation_matrix.bottomLeftCorner(1,3) = center->transpose();

  algebraic_form = translation_matrix * algebraic_form * translation_matrix.transpose();

  //std::cout << "translation matrix: \n" << translation_matrix << std::endl;
  SelfAdjointEigenSolver< Matrix< ScalarType, 3, 3 > > eigen_solver( algebraic_form.topLeftCorner(3,3)
      * (1.0/-algebraic_form(3,3)) );

  *eigenvalue_matrix = eigen_solver.eigenvalues().asDiagonal();
  std::cout << "eigen values: \n" << *eigenvalue_matrix << std::endl;

  *eigenvectors = eigen_solver.eigenvectors();
  std::cout << "eigen vectors: \n" << *eigenvectors << std::endl;

}

//generate random unit sphere points and distort them according to an fixed ellipsoid
//rotation, translation and scaling matrix and some white noise
template< typename ScalarType, unsigned int number_of_points>
void GenerateTestValues( Array< ScalarType, 3, number_of_points> *result_sphere,
    Array< ScalarType, 3, number_of_points> *result_ellipsoid,
    const ScalarType white_noise_variance )
{

  //model soft iron, axis scaling error
  //transform and scale points

  //scale error
  Matrix< ScalarType, 3, 3> S_M;
  S_M <<  1.2, 0.0, 0.0,
          0.0, 0.8, 0.0,
          0.0, 0.0, 1.3;

  //soft iron error
  Matrix< ScalarType, 3, 3> C_SI;
  C_SI << 0.58, -0.73, 0.36,
          1.32, 0.46, -0.12,
          -0.26, 0.44, 0.53;

  //hard iron offset
  Matrix< ScalarType, 3, 1> b_HI;
  b_HI << -1.2,
           0.2,
          -0.8;

  //axis offset
  Matrix< ScalarType, 3, 1> b_M;
  b_M <<  1.5,
          0.4,
          2.7;
  

  std::default_random_engine generator;
  for (unsigned int i = 0; i < number_of_points; ++i)
  {
    Matrix< ScalarType, 3, 1> random_sphere_point;

    std::uniform_real_distribution< ScalarType > distribution_x(-1.0, 1.0);

    //generate unit sphere points according to equation
    //x²+y²+z² = 1;

    //generate x coordinate
    random_sphere_point.x() = distribution_x(generator);

    //generate y coordinate
    //rescale distribution values according the random_sphere_point.x
    ScalarType a = -(1 - (random_sphere_point.x()*random_sphere_point.x()));
    ScalarType b = 1 - random_sphere_point.x()*random_sphere_point.x();

    std::uniform_real_distribution< ScalarType > distribution_y(a, b);
    random_sphere_point.y() = distribution_y(generator);

    //generate z coordinate
    a = distribution_y.a() + (random_sphere_point.y()*random_sphere_point.y());
    b = distribution_y.b() - (random_sphere_point.y()*random_sphere_point.y());

    std::uniform_real_distribution< ScalarType > distribution_z(a, b);
    random_sphere_point.z() = distribution_z(generator);

    //scale vector to unit sphere point
    random_sphere_point *= 1.0/random_sphere_point.norm();

    result_sphere->col(i) = random_sphere_point;

    //distort sphere point
    result_ellipsoid->col(i) = C_SI*S_M * random_sphere_point + b_HI + b_M;

    //TODO add noise

    //std::cout << "Vector: " << result_ellipsoid->col(i) << std::endl; 
    //std::cout << "Norm: " << random_sphere_point.norm() << std::endl; 
  }
  std::cout << "real center vector: " << b_HI + b_M << std::endl; 

}


