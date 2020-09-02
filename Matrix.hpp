//
// Created by amirt01 on 9/1/20.
//

#ifndef CPP_MATRIX_MANIPULATION_MATRIX_HPP
#define CPP_MATRIX_MANIPULATION_MATRIX_HPP

#include <ostream>
#include <vector>
#include <algorithm>

template<class T>
class Matrix {
 private:
  std::vector<std::vector<T>> m_matrix;
  unsigned m_rows = 0;
  unsigned m_cols = 0;
 public:
  Matrix() = default;

  Matrix(const unsigned& rows, const unsigned& cols, const T& initial = T{});

  explicit Matrix(const std::vector<std::vector<T>>& matrix);

  Matrix(const std::initializer_list<std::vector<T>>& il);

  Matrix(const Matrix<T>& rhs);

  void resize(const unsigned& rows, const unsigned& cols, const T& initial = T{});

  // Operator overloading, for "standard" mathematical matrix operations
  Matrix<T>& operator=(const Matrix<T>& rhs);

  // Matrix mathematical operations
  Matrix<T> operator+(const Matrix<T>& rhs);

  Matrix<T>& operator+=(const Matrix<T>& rhs);

  Matrix<T> operator-(const Matrix<T>& rhs);

  Matrix<T>& operator-=(const Matrix<T>& rhs);

  Matrix<T> operator*(const Matrix<T>& rhs);

  Matrix<T>& operator*=(const Matrix<T>& rhs);

  void transpose();

  // Matrix/scalar operations
  Matrix<T> operator+(const T& rhs);

  Matrix<T> operator-(const T& rhs);

  Matrix<T> operator*(const T& rhs);

  Matrix<T> operator/(const T& rhs);

  Matrix<T> operator+=(const T& rhs);

  Matrix<T> operator-=(const T& rhs);

  Matrix<T> operator*=(const T& rhs);

  Matrix<T> operator/=(const T& rhs);

  // Matrix/vector operations
  std::vector<T> operator*(const std::vector<T>& rhs);

  std::vector<T> diag_vec();

  // Access the individual elements
  T& operator()(const unsigned& row, const unsigned& col);

  const T& operator()(const unsigned& row, const unsigned& col) const;

  // Access the row and column sizes
  [[nodiscard]] unsigned get_rows() const;

  [[nodiscard]] unsigned get_cols() const;

  // Insertion operator overload
  friend std::ostream& operator<<(std::ostream& os, const Matrix<T> rhs) {
    // print each cell in order
    for (const auto& row : rhs.m_matrix) {
      for (const auto& cell : row) {
        os << cell << ' ';
      }
      // separate each row
      os << std::endl;
    }
    return os;
  }
};

// Parameter Constructor
template<class T>
Matrix<T>::Matrix(const unsigned& rows, const unsigned& cols, const T& initial)
  : m_rows(rows), m_cols(cols) {
  resize(rows, cols, initial);
}

// Matrix initializer
template<class T>
Matrix<T>::Matrix(const std::vector<std::vector<T>>& matrix)
  : m_matrix(matrix), m_rows(matrix.size()), m_cols(matrix[0].size()) {}

// Vector initializer list
template<class T>
Matrix<T>::Matrix(const std::initializer_list<std::vector<T>>& il)
  : m_matrix(il), m_rows(il.size()), m_cols(il.begin()->size()) {
  // check to make sure the matrix is all the correct size
  for (const auto& row : m_matrix) {
    if (row.size() != m_cols) {
      throw "INVALID MATRIX";
    }
  }
}

// Copy Constructor
template<class T>
Matrix<T>::Matrix(const Matrix<T>& rhs)
  : m_matrix(rhs.m_matrix), m_rows(rhs.m_rows), m_cols(rhs.m_cols) {}

template<class T>
void Matrix<T>::resize(const unsigned int& rows, const unsigned int& cols, const T& initial) {
  // resize the number of rows in the matrix
  m_matrix.resize(rows);
  // resize the number of columns in the matrix
  for (auto& col : m_matrix) {
    col.resize(cols, initial);
  }
}

// Assignment Operator
template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs) {
  // check if the matrices are already equal
  if (&rhs == this)
    return *this;

  // copy the rows and cols from rhs
  m_rows = rhs.m_rows;
  m_cols = rhs.m_cols;

  // copy the matrix from rhs
  m_matrix = rhs.m_matrix;

  return *this;
}

// Addition of two matrices
template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs) {
  // test if matrices are the same size
  if (rhs.m_rows != m_rows) {
    throw "MISMATCHING NUMBER OF ROWS";
  }
  if (rhs.m_cols != m_cols) {
    throw "MISMATCHING NUMBER OF COLS";
  }

  Matrix result(m_rows, m_cols);

  // add each cell from this matrix and rhs into the resulting matrix
  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      result(i, j) = m_matrix[i][j] + rhs(i, j);
    }
  }

  // move ownership in function return
  return std::move(result);
}

// Cumulative addition of this matrix and another
template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {
  // test if matrices are the same size
  if (rhs.m_rows != m_rows) {
    throw "MISMATCHING NUMBER OF ROWS";
  }
  if (rhs.m_cols != m_cols) {
    throw "MISMATCHING NUMBER OF COLS";
  }

  // add each cell from rhs to matrix
  for (unsigned i = 0; i < rhs.m_rows; i++) {
    for (unsigned j = 0; j < rhs.m_cols; j++) {
      this->m_matrix[i][j] += rhs(i, j);
    }
  }

  return *this;
}

// Subtraction of this matrix and another
template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) {
  // test if matrices are the same size
  if (rhs.m_rows != m_rows) {
    throw "MISMATCHING NUMBER OF ROWS";
  }
  if (rhs.m_cols != m_cols) {
    throw "MISMATCHING NUMBER OF COLS";
  }

  Matrix result(m_rows, m_cols);

  // subtract each cell of rhs from this matrix and store in the resulting matrix
  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      result(i, j) = this->m_matrix[i][j] - rhs(i, j);
    }
  }

  // move ownership in function return
  return std::move(result);
}

// Cumulative subtraction of this matrix and another
template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {
  // test if matrices are the same size
  if (rhs.m_rows != m_rows) {
    throw "MISMATCHING NUMBER OF ROWS";
  }
  if (rhs.m_cols != m_cols) {
    throw "MISMATCHING NUMBER OF COLS";
  }

  // subtract each cell of rhs from matrix
  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      this->mat[i][j] -= rhs(i, j);
    }
  }

  return *this;
}

// Left multiplication of this matrix and another
template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) {
  // test if rows of the second matrix is equal to the columns of the other
  if (rhs.m_rows != m_cols) {
    throw "MISMATCHING NUMBER OF ROWS";
  }

  Matrix result(m_rows, m_cols);

  // multiply the matrices
  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      for (unsigned k = 0; k < m_rows; k++) {
        result(i, j) += this->m_matrix[i][k] * rhs(k, j);
      }
    }
  }

  return std::move(result);
}

// Cumulative left multiplication of this matrix and another
template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs) {
  Matrix result = *this * rhs;
  *this = result;
  return *this;
}

// Transpose this matrix (square)
template<class T>
void Matrix<T>::transpose() {
  // check if the matrix is square
  if (m_cols != m_rows) {
    throw "MATRIX NOT SQUARE";
  }

  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < i; j++) {
      std::swap(m_matrix[i][j], m_matrix[j][i]);
    }
  }
}

// Matrix/scalar addition
template<class T>
Matrix<T> Matrix<T>::operator+(const T& rhs) {
  Matrix result(m_rows, m_cols);

  // add rhs to each cell in the matrix
  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      result(i, j) = this->m_matrix[i][j] + rhs;
    }
  }

  return std::move(result);
}

// Matrix/scalar subtraction
template<class T>
Matrix<T> Matrix<T>::operator-(const T& rhs) {
  Matrix result(m_rows, m_cols);

  // subtract rhs from each cell in the matrix
  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      result(i, j) = this->m_matrix[i][j] - rhs;
    }
  }

  return std::move(result);
}

// Matrix/scalar multiplication
template<typename T>
Matrix<T> Matrix<T>::operator*(const T& rhs) {
  Matrix result(m_rows, m_cols);

  // multiply each cell in the matrix by rhs
  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      result(i, j) = this->m_matrix[i][j] * rhs;
    }
  }

  return std::move(result);
}

// Matrix/scalar division
template<typename T>
Matrix<T> Matrix<T>::operator/(const T& rhs) {
  Matrix result(m_rows, m_cols, 0.0);

  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      result(i, j) = this->m_matrix[i][j] / rhs;
    }
  }

  return std::move(result);
}

template<class T>
Matrix<T> Matrix<T>::operator+=(const T& rhs) {
  for (auto& row : m_matrix) {
    for (auto& cell : row) {
      cell += rhs;
    }
  }
}

template<class T>
Matrix<T> Matrix<T>::operator-=(const T& rhs) {
  for (auto& row : m_matrix) {
    for (auto& cell : row) {
      cell -= rhs;
    }
  }
}

template<class T>
Matrix<T> Matrix<T>::operator*=(const T& rhs) {
  for (auto& row : m_matrix) {
    for (auto& cell : row) {
      cell *= rhs;
    }
  }
}

template<class T>
Matrix<T> Matrix<T>::operator/=(const T& rhs) {
  for (auto& row : m_matrix) {
    for (auto& cell : row) {
      cell /= rhs;
    }
  }
}

// Multiply a matrix with a vector
template<class T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& rhs) {
  std::vector<T> result(rhs.size(), 0.0);

  for (unsigned i = 0; i < m_rows; i++) {
    for (unsigned j = 0; j < m_cols; j++) {
      result[i] = this->m_matrix[i][j] * rhs[j];
    }
  }

  return std::move(result);
}

// Obtain a vector of the diagonal elements
template<class T>
std::vector<T> Matrix<T>::diag_vec() {
  std::vector<T> result(m_rows, 0.0);

  for (unsigned i = 0; i < m_rows; i++) {
    result[i] = this->m_matrix[i][i];
  }

  return std::move(result);
}

// Access the individual elements
template<class T>
T& Matrix<T>::operator()(const unsigned int& row, const unsigned int& col) {
  return m_matrix[row][col];
}

// Access the individual elements (const)
template<typename T>
const T& Matrix<T>::operator()(const unsigned& row, const unsigned& col) const {
  return m_matrix[row][col];
}

// Get the number of rows of the matrix
template<typename T>
unsigned Matrix<T>::get_rows() const {
  return m_rows;
}

// Get the number of columns of the matrix
template<typename T>
unsigned Matrix<T>::get_cols() const {
  return m_cols;
}

#endif //CPP_MATRIX_MANIPULATION_MATRIX_HPP
