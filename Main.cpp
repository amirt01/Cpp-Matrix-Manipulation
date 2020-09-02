//
// Created by amirt01 on 9/1/20.
//

#include <iostream>
#include "Matrix.hpp"

int main() {
  Matrix<double> mat{{1, 2, 3},
                     {4, 5, 6}};

  std::cout << mat;

  std::cout << '\n';
  mat.transpose();

  std::cout << mat;


  std::cout << '\n';
  mat.resize(3, 4, 5);

  std::cout << mat;

  return EXIT_SUCCESS;
}
