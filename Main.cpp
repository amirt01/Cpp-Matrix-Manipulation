//
// Created by amirt01 on 9/1/20.
//

#include <iostream>
#include "Matrix.hpp"

int main() {
  Matrix<double> mat{{1, 2, 3},
                     {4, 5, 6}};

  std::cout << "ORIGINAL VARIABLE:\n";
  std::cout << mat;

  std::cout << '\n';

  std::cout << "TRANSPOSED MATRIX:\n";
  std::cout << transpose(mat);

  std::cout << '\n';

  std::cout << "ORIGINAL VARIABLE:\n";
  std::cout << mat;

  return EXIT_SUCCESS;
}
