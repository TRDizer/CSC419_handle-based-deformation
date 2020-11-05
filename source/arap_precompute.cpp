#include "../include/arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/arap_linear_block.h>
#include <igl/cotmatrix.h>
#include <iostream>

#define get_ops_e1(e)    ((e + 1) % 3)
#define get_ops_e2(e)    ((e + 2) % 3)

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
  // REPLACE WITH YOUR CODE
  Eigen::SparseMatrix<double> lapacian;
  igl::cotmatrix(V, F, lapacian);

  igl::min_quad_with_fixed_precompute(lapacian, b, Eigen::SparseMatrix<double>(), false, data);
  
  K.resize(V.rows(), 3 * V.rows());
  K.setZero();
  std::vector<Eigen::Triplet<double>> K_entries;

  int ei_index, ej_index, k_entry;
  for (int face_index = 0; face_index < F.rows(); face_index++) {
      for (int edge = 0; edge < 3; edge++) {
          // loop over to get each edge difference 
          ei_index = F(face_index, get_ops_e1(edge));
          ej_index = F(face_index, get_ops_e2(edge));

          Eigen::RowVector3d tilde_e_diff = V.row(ei_index) - V.row(ej_index);
          tilde_e_diff *= lapacian.coeff(ei_index, ej_index) / 6.;

          for (int vertex = 0; vertex < 3; vertex++) {
              // loop over to populate per-vertex rotation
              k_entry = F(face_index, vertex);

              // Seems like there is no block write for Eigen sparse matrix
              // K.block<1,3>(ei_index, 3 * k) = tilde_e_diff
              // K.block<1,3>(ej_index, 3 * k) = tilde_e_diff

              for (int beta = 0; beta < 3; beta++) {
                // K_i_(3k + beta) += tilde_eij^beta 
                // K_j_(3k + beta) -= tilde_eij^beta 
                K_entries.push_back(Eigen::Triplet<double>(ei_index, 3 * k_entry + beta, tilde_e_diff(beta)));
                K_entries.push_back(Eigen::Triplet<double>(ej_index, 3 * k_entry + beta, -tilde_e_diff(beta)));
              }
          }
      }
  }

  K.setFromTriplets(K_entries.begin(), K_entries.end());
  std::cout << "# non-zero entries is " << K.nonZeros() << std::endl;
}
