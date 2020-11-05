#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
  // REPLACE WITH YOUR CODE
  // data.n = V.rows();

  Eigen::SparseMatrix<double> lapacian;
  igl::cotmatrix(V, F, lapacian);

  Eigen::SparseMatrix<double> mass, mass_inv, Q;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);
  igl::invert_diag(mass, mass_inv);
  Q = lapacian.transpose() * mass_inv * lapacian;

  igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);
}
