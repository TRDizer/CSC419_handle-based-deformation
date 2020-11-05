#include "../include/arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
  // REPLACE WITH YOUR CODE

  // local step tr(V^T * K * R) where V^T * K gives C^T
  Eigen::MatrixXd C = (K.transpose()*U).normalized();
  Eigen::MatrixXd R(3 * U.rows(), 3);

  Eigen::Matrix3d Rk;
  for(int index = 0; index < U.rows(); index++){
    // Solve all diagonal blocks to obtain the closest rotation
    igl::polar_svd3x3(Eigen::Matrix3d(C.block<3,3>(3 * index,0)), Rk);
    R.block<3,3>(3 * index, 0) = Rk;
  }

  //global step: quadratic minimization on displacement term since we have rotation now
  Eigen::MatrixXd B = K*R;
  igl::min_quad_with_fixed_solve(data, B, bc, Eigen::VectorXd(), U);
}
