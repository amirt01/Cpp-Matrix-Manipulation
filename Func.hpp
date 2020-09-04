//
// Created by amirt01 on 9/2/20.
//

#ifndef CPP_MATRIX_MANIPULATION_FUNC_HPP
#define CPP_MATRIX_MANIPULATION_FUNC_HPP

#include "Matrix.hpp"

// transpose the matrix without
template<typename T>
Matrix<T> transpose(Matrix<T> rhs) {
  Matrix<T> result(rhs.get_cols(), rhs.get_rows());
  for (unsigned i = 0; i < rhs.get_rows(); i++) {
    for (unsigned j = 0; j < rhs.get_cols(); j++) {
      result(j, i) = rhs(i, j);
    }
  }
  return std::move(result);
}

#endif //CPP_MATRIX_MANIPULATION_FUNC_HPP
