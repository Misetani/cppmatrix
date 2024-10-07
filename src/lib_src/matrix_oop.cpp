#include <cmath>
#include <cstdio>

#include "../array_exception.h"
#include "../s21_matrix_oop.h"

/* ========================================================================= */
/*                       Constructors and Destructors                        */
/* ========================================================================= */

S21Matrix::S21Matrix(int r, int c) : rows{r}, cols{c} {
  matrix = new double*[rows];

  for (int i = 0; i < rows; ++i) {
    matrix[i] = new double[cols]{0};
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows, other.cols) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      matrix[i][j] = other.matrix[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows{other.rows}, cols{other.cols}, matrix{other.matrix} {
  // Steal the resources from 'other'
  other.matrix = nullptr;
  other.rows = 0;
  other.cols = 0;
}

S21Matrix::~S21Matrix() {
  free();

  rows = 0;
  cols = 0;
}

/* ========================================================================= */
/*                            Getters and Setters                            */
/* ========================================================================= */

// trying to be efficient gives me a lot of pain
void S21Matrix::set_rows(int r) {
  if (r < 0) {
    throw Array_exception("Number of rows must be positive");
  }

  if (r < rows) {  // free excessive rows
    double** m = new double*[r];

    for (int i = 0; i < r; ++i) {
      m[i] = matrix[i];
    }

    for (int i = r; i < rows; ++i) {
      delete[] matrix[i];
    }

    delete[] matrix;

    matrix = m;
  } else if (r > rows) {  // allocate enough memory for new rows
    double** m = new double*[r];

    for (int i = 0; i < rows; ++i) {
      m[i] = matrix[i];
    }

    for (int i = rows; i < r; ++i) {
      m[i] = new double[cols]{0};
    }

    delete[] matrix;

    matrix = m;
  }

  rows = r;
}

void S21Matrix::set_cols(int c) {
  if (c < 0) {
    throw Array_exception("Number of columns must be positive!");
  }

  double** m = new double*[rows];

  for (int i = 0; i < rows; ++i) {
    m[i] = new double[c]{0};

    for (int j = 0; j < (c < cols ? c : cols); ++j) {
      m[i][j] = matrix[i][j];
    }
  }

  free();

  cols = c;

  matrix = m;
}

/* ========================================================================= */
/*                               Main Methods                                */
/* ========================================================================= */

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (matrix == nullptr || other.matrix == nullptr) {
    return false;
  }

  if (rows != other.rows || cols != other.cols) {
    return false;
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (std::fabs(matrix[i][j] - other.matrix[i][j]) > 1e-7) {
        return false;
      }
    }
  }

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows != other.rows || cols != other.cols) {
    throw Array_exception("Cannot add matrices with different dimensions!");
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      matrix[i][j] += other.matrix[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows != other.rows || cols != other.cols) {
    throw Array_exception(
        "Cannot substract from matrix with different dimensions!");
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      matrix[i][j] -= other.matrix[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      matrix[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols != other.rows) {
    throw Array_exception(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix");
  }

  // product of two matrixces is of different dimensions than
  // any of the matrices being multiplied in all cases except
  // for the square matrices
  S21Matrix first = *this;

  init_matrix(first.rows, other.cols);

  for (int row = 0; row < first.rows; ++row) {
    for (int column = 0; column < other.cols; ++column) {
      for (int i = 0; i < first.cols; ++i) {
        matrix[row][column] += first.matrix[row][i] * other.matrix[i][column];
      }
    }
  }
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix result(cols, rows);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result.matrix[j][i] = matrix[i][j];
    }
  }

  return result;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows != cols) {
    throw Array_exception("Matrix is not a square matrix!");
  }

  S21Matrix result(cols, rows);

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      S21Matrix sub_matrix(rows - 1, cols - 1);

      for (int k = 0; k < rows; ++k) {
        if (k < i) {
          for (int l = 0; l < cols; ++l) {
            if (l < j) {
              sub_matrix.matrix[k][l] = matrix[k][l];
            } else if (l > j) {
              sub_matrix.matrix[k][l - 1] = matrix[k][l];
            }
          }
        } else if (k > i) {
          for (int l = 0; l < cols; ++l) {
            if (l < j) {
              sub_matrix.matrix[k - 1][l] = matrix[k][l];
            } else if (l > j) {
              sub_matrix.matrix[k - 1][l - 1] = matrix[k][l];
            }
          }
        }
      }

      result.matrix[i][j] = sub_matrix.Determinant();

      result.matrix[i][j] *= (i + j) % 2 == 0 ? 1 : -1;
    }
  }

  return result;
}

double S21Matrix::Determinant() const {
  if (rows != cols) {
    throw Array_exception("Matrix is not a square matrix!");
  }

  if (rows == 1) {  // base case
    return matrix[0][0];
  }

  if (rows == 2) {  // base case
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  }

  double result = 0.0;

  // iterate over the first row
  for (int i = 0; i < cols; ++i) {
    S21Matrix sub_matrix(rows - 1, cols - 1);

    for (int j = 1; j < rows; ++j) {
      for (int k = 0; k < cols; ++k) {
        if (k < i) {
          sub_matrix.matrix[j - 1][k] = matrix[j][k];
        } else if (k > i) {
          sub_matrix.matrix[j - 1][k - 1] = matrix[j][k];
        }
      }
    }

    int sign = (i % 2 == 0) ? 1 : -1;

    double sub_det = sub_matrix.Determinant();

    result += sign * matrix[0][i] * sub_det;
  }

  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (rows != cols) {
    throw Array_exception(
        "Cannot calculate determinant for a non-squre matrix!");
  }

  double det = Determinant();

  if (det == 0.0) {
    throw Array_exception("Matrix determinant is 0!");
  }

  S21Matrix result(cols, rows);

  if (rows == 1) {
    result.matrix[0][0] = 1 / matrix[0][0];
    return result;
  }

  S21Matrix transposed_complements = CalcComplements().Transpose();

  for (int i = 0; i < result.rows; ++i) {
    for (int j = 0; j < result.cols; ++j) {
      result.matrix[i][j] = transposed_complements.matrix[i][j] / det;
    }
  }

  return result;
}

/* ========================================================================= */
/*                             Matrix Operators                              */
/* ========================================================================= */

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix result = *this;

  result.SumMatrix(other);

  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix result = *this;

  result.SubMatrix(other);

  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix result = *this;

  result.MulMatrix(other);

  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return this->EqMatrix(other);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  init_matrix(other.rows, other.cols);

  for (int i = 0; i < other.rows; ++i) {
    for (int j = 0; j < other.rows; ++j) {
      matrix[i][j] = other.matrix[i][j];
    }
  }

  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(double num) {
  this->MulNumber(num);
  return *this;
}

/* ========================================================================= */
/*                          Helper public methods                            */
/* ========================================================================= */

double& S21Matrix::get_element(int row, int col) {
  if (row < 0 || row >= rows) {
    throw Array_exception("Row indexation is out of bonds!");
  }

  if (col < 0 || col >= cols) {
    throw Array_exception("Column indexation is out of bonds!");
  }

  return matrix[row][col];
}

const double& S21Matrix::get_element(int row, int col) const {
  if (row < 0 || row >= rows) {
    throw Array_exception("Row indexation is out of bonds!");
  }

  if (col < 0 || col >= cols) {
    throw Array_exception("Column indexation is out of bonds!");
  }

  return matrix[row][col];
}

// void S21Matrix::print() const {
//     for (int i = 0; i < rows; ++i) {
//         for (int j = 0; j < cols; ++j) {
//             printf("%-7.2lf", get_element(i, j));
//         }
//         printf("\n");
//     }
// }

/* ========================================================================= */
/*                          Helper private methods                           */
/* ========================================================================= */

void S21Matrix::free() {
  if (matrix == nullptr) {
    return;
  }

  for (int i = 0; i < rows; ++i) {
    delete[] matrix[i];  // delete[] frees the memory allocated to a block
                         // starting with specified pointer
  }
  delete[] matrix;

  matrix = nullptr;
}

void S21Matrix::init_matrix(int r, int c) {
  free();

  rows = r;
  cols = c;

  matrix = new double*[rows];

  for (int i = 0; i < rows; ++i) {
    matrix[i] = new double[cols]{0};
  }
}