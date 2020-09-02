//
// Created by amirt01 on 9/1/20.
//

#include <iostream>
#include "Matrix.hpp"
#include "Func.hpp"

int main() {
  Matrix<double> mat{{1, 2, 3},
                     {4, 5, 6}};

  std::cout << mat;

  std::cout << '\n';

  std::cout << transpose(mat);

  std::cout << 'n';

  std::cout << mat;

  return EXIT_SUCCESS;
}
