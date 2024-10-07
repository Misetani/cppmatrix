#include <gtest/gtest.h>

#include "../array_exception.h"
#include "../s21_matrix_oop.h"

S21Matrix get_basic_matrix(int rows, int cols) {
  S21Matrix matrix(rows, cols);

  int count = 1;
  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      matrix.get_element(i, j) = count++;
    }
  }

  return matrix;
}

S21Matrix get_basic_sum(int rows, int cols) {
  S21Matrix matrix(rows, cols);

  int count = 1;
  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      matrix.get_element(i, j) = 2 * count;
      count += 1;
    }
  }

  return matrix;
}

S21Matrix get_basic_dif(int rows, int cols) {
  S21Matrix matrix(rows, cols);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      matrix.get_element(i, j) = 0;
    }
  }

  return matrix;
}

S21Matrix get_matrix_32() {
  S21Matrix matrix(3, 2);

  matrix.get_element(0, 0) = 1.0;
  matrix.get_element(0, 1) = 2.0;
  matrix.get_element(1, 0) = 3.0;
  matrix.get_element(1, 1) = 4.0;
  matrix.get_element(2, 0) = 5.0;
  matrix.get_element(2, 1) = 6.0;

  return matrix;
}

S21Matrix get_matrix_23() {
  S21Matrix matrix(2, 3);

  matrix.get_element(0, 0) = 1.0;
  matrix.get_element(0, 1) = 2.0;
  matrix.get_element(0, 2) = 3.0;
  matrix.get_element(1, 0) = 4.0;
  matrix.get_element(1, 1) = 5.0;
  matrix.get_element(1, 2) = 6.0;

  return matrix;
}

S21Matrix get_basic_product() {
  S21Matrix matrix(3, 3);

  matrix.get_element(0, 0) = 9.0;
  matrix.get_element(0, 1) = 12.0;
  matrix.get_element(0, 2) = 15.0;
  matrix.get_element(1, 0) = 19.0;
  matrix.get_element(1, 1) = 26.0;
  matrix.get_element(1, 2) = 33.0;
  matrix.get_element(2, 0) = 29.0;
  matrix.get_element(2, 1) = 40.0;
  matrix.get_element(2, 2) = 51.0;

  return matrix;
}

S21Matrix get_matrix_for_complements() {
  S21Matrix A(3, 3);

  A.get_element(0, 0) = 1.0;

  A.get_element(0, 0) = 1.0;
  A.get_element(0, 1) = 2.0;
  A.get_element(0, 2) = 3.0;
  A.get_element(1, 0) = 0.0;
  A.get_element(1, 1) = 4.0;
  A.get_element(1, 2) = 2.0;
  A.get_element(2, 0) = 5.0;
  A.get_element(2, 1) = 2.0;
  A.get_element(2, 2) = 1.0;

  return A;
}

S21Matrix get_matrix_of_compliments() {
  S21Matrix A(3, 3);

  A.get_element(0, 0) = 0.0;
  A.get_element(0, 1) = 10.0;
  A.get_element(0, 2) = -20.0;
  A.get_element(1, 0) = 4.0;
  A.get_element(1, 1) = -14.0;
  A.get_element(1, 2) = 8.0;
  A.get_element(2, 0) = -8.0;
  A.get_element(2, 1) = -2.0;
  A.get_element(2, 2) = 4.0;

  return A;
}

S21Matrix get_matrix_for_inversion() {
  S21Matrix A(3, 3);

  A.get_element(0, 0) = 2.0;
  A.get_element(0, 1) = 5.0;
  A.get_element(0, 2) = 7.0;
  A.get_element(1, 0) = 6.0;
  A.get_element(1, 1) = 3.0;
  A.get_element(1, 2) = 4.0;
  A.get_element(2, 0) = 5.0;
  A.get_element(2, 1) = -2.0;
  A.get_element(2, 2) = -3.0;

  return A;
}

S21Matrix get_inverted_matrix() {
  S21Matrix A(3, 3);

  A.get_element(0, 0) = 1.0;
  A.get_element(0, 1) = -1.0;
  A.get_element(0, 2) = 1.0;
  A.get_element(1, 0) = -38.0;
  A.get_element(1, 1) = 41.0;
  A.get_element(1, 2) = -34.0;
  A.get_element(2, 0) = 27.0;
  A.get_element(2, 1) = -29.0;
  A.get_element(2, 2) = 24.0;

  return A;
}

TEST(creating_a_matrix, default_matrix) {
  S21Matrix matrix;  // default constructor should be called;
  EXPECT_EQ(matrix.get_rows(), 0);
  EXPECT_EQ(matrix.get_cols(), 0);
}

TEST(creating_a_matrix, simple_matrix) {
  S21Matrix matrix(3, 3);
  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 3);

  int count = 1;
  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      matrix.get_element(i, j) = count++;
    }
  }

  count = 1;
  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), count++);
    }
  }

  EXPECT_THROW(matrix.get_element(-1, 0), Array_exception);
  EXPECT_THROW(matrix.get_element(0, -1), Array_exception);
  EXPECT_THROW(matrix.get_element(3, 0), Array_exception);
  EXPECT_THROW(matrix.get_element(0, 3), Array_exception);
  EXPECT_THROW(matrix.get_element(999, 999), Array_exception);
}

TEST(changing_matrix_dimensions, set_rows) {
  S21Matrix matrix(3, 3);
  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 3);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), 0);
    }
  }

  matrix.set_rows(7);

  EXPECT_EQ(matrix.get_rows(), 7);
  EXPECT_EQ(matrix.get_cols(), 3);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), 0);
    }
  }

  matrix.set_rows(2);

  EXPECT_EQ(matrix.get_rows(), 2);
  EXPECT_EQ(matrix.get_cols(), 3);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), 0);
    }
  }

  EXPECT_THROW(matrix.set_rows(-1), Array_exception);
}

TEST(changing_matrix_dimensions, set_cols) {
  S21Matrix matrix(3, 3);
  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 3);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), 0);
    }
  }

  matrix.set_cols(7);

  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 7);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), 0);
    }
  }

  matrix.set_cols(2);

  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 2);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), 0);
    }
  }

  EXPECT_THROW(matrix.set_cols(-1), Array_exception);
  EXPECT_THROW(matrix.set_rows(-1), Array_exception);
}

TEST(copying_a_matrix, with_copy_constructor) {
  S21Matrix matrix(3, 3);
  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 3);

  int count = 1;
  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      matrix.get_element(i, j) = count++;
    }
  }

  S21Matrix matrix2 = matrix;

  EXPECT_EQ(matrix.get_rows(), matrix2.get_rows());
  EXPECT_EQ(matrix.get_cols(), matrix2.get_cols());

  count = 1;
  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix2.get_element(i, j), matrix.get_element(i, j));
    }
  }
}

TEST(copying_a_matrix, with_move_constructor) {
  // calling a move constructor is not so easy because of an ample optimisation
  S21Matrix matrix0 = S21Matrix(3, 3);
  S21Matrix matrix = std::move(matrix0);

  EXPECT_EQ(matrix.get_rows(), 3);
  EXPECT_EQ(matrix.get_cols(), 3);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j), 0);
    }
  }
}

TEST(comparing_matrices, on_equality) {
  S21Matrix A(3, 3);
  S21Matrix B(3, 3);

  EXPECT_TRUE(A.EqMatrix(B));

  A.get_element(1, 1) = 1.0;
  EXPECT_FALSE(A.EqMatrix(B));

  B.get_element(1, 1) = 1.000001;
  EXPECT_FALSE(A.EqMatrix(B));

  B.get_element(1, 1) = 1.0000001;
  EXPECT_FALSE(A.EqMatrix(B));

  B.get_element(1, 1) = 1.00000001;
  EXPECT_TRUE(A.EqMatrix(B));
}

TEST(comparing_matrices, on_equality_2) {
  S21Matrix A;
  S21Matrix B;

  EXPECT_FALSE(A.EqMatrix(B));

  S21Matrix C(2, 3);
  S21Matrix D(3, 2);

  EXPECT_FALSE(A.EqMatrix(C));
  EXPECT_FALSE(A.EqMatrix(D));
  EXPECT_FALSE(C.EqMatrix(D));
}

TEST(matrix_arithmetic, matrix_addition) {
  S21Matrix matrix = get_basic_matrix(3, 3);

  matrix.SumMatrix(get_basic_matrix(3, 3));

  EXPECT_TRUE(matrix.EqMatrix(get_basic_sum(3, 3)));
}

TEST(matrix_arithmetic, matrix_addition_exceptions) {
  S21Matrix A = get_basic_matrix(2, 3);
  S21Matrix B = get_basic_matrix(3, 3);
  S21Matrix C = get_basic_matrix(3, 2);

  EXPECT_THROW(A.SumMatrix(B), Array_exception);
  EXPECT_THROW(A.SumMatrix(C), Array_exception);
  EXPECT_THROW(B.SumMatrix(C), Array_exception);
}

TEST(matrix_arithmetic, matrix_subtraction) {
  S21Matrix matrix = get_basic_matrix(3, 3);

  matrix.SubMatrix(get_basic_matrix(3, 3));

  EXPECT_TRUE(matrix.EqMatrix(get_basic_dif(3, 3)));
}

TEST(matrix_arithmetic, matrix_subtraction_exceptions) {
  S21Matrix A = get_basic_matrix(2, 3);
  S21Matrix B = get_basic_matrix(3, 3);
  S21Matrix C = get_basic_matrix(3, 2);

  EXPECT_THROW(A.SubMatrix(B), Array_exception);
  EXPECT_THROW(A.SubMatrix(C), Array_exception);
  EXPECT_THROW(B.SubMatrix(C), Array_exception);
}

TEST(matrix_arithmetic, matrix_number_multiplication) {
  S21Matrix matrix = get_basic_matrix(3, 3);
  S21Matrix basic_matrix = get_basic_matrix(3, 3);

  matrix.MulNumber(3.14);

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j),
                basic_matrix.get_element(i, j) * 3.14);
    }
  }
}

TEST(matrix_arithmetic, matrix_multiplication) {
  S21Matrix matrix1 = get_matrix_32();
  S21Matrix matrix2 = get_matrix_23();

  matrix1.MulMatrix(matrix2);

  EXPECT_TRUE(matrix1.EqMatrix(get_basic_product()));
}

TEST(matrix_arithmetic, matrix_multiplication_exceptions) {
  S21Matrix matrix1 = get_matrix_23();
  S21Matrix matrix2 = get_matrix_23();

  EXPECT_THROW(matrix2.MulMatrix(matrix1), Array_exception);
}

TEST(matrix_transformation, matrix_determinant_1) {
  S21Matrix matrix(1, 1);
  matrix.get_element(0, 0) = 1.0;

  EXPECT_NEAR(matrix.Determinant(), 1.0, 1e-7);
}

TEST(matrix_transformation, matrix_determinant_2) {
  S21Matrix matrix(2, 2);
  matrix.get_element(0, 0) = 1.0;
  matrix.get_element(0, 1) = 2.0;
  matrix.get_element(1, 0) = 3.0;
  matrix.get_element(1, 1) = 4.0;

  EXPECT_NEAR(matrix.Determinant(), -2.0, 1e-7);
}

TEST(matrix_transformation, matrix_determinant_3) {
  S21Matrix matrix(3, 3);
  matrix.get_element(0, 0) = 1.0;
  matrix.get_element(0, 1) = 2.0;
  matrix.get_element(0, 2) = 3.0;
  matrix.get_element(1, 0) = 4.0;
  matrix.get_element(1, 1) = 5.0;
  matrix.get_element(1, 2) = 6.0;
  matrix.get_element(2, 0) = 7.0;
  matrix.get_element(2, 1) = 8.0;
  matrix.get_element(2, 2) = 9.0;

  EXPECT_NEAR(matrix.Determinant(), 0.0, 1e-7);
}

TEST(matrix_transformation, matrix_determinant_4) {
  S21Matrix matrix(3, 4);

  EXPECT_THROW(matrix.Determinant(), Array_exception);
}

TEST(matrix_transformation, matrix_transposition) {
  S21Matrix matrix = get_basic_matrix(3, 3);

  matrix = matrix.Transpose();

  int count = 1;
  for (int j = 0; j < matrix.get_cols(); ++j) {
    for (int i = 0; i < matrix.get_rows(); ++i) {
      EXPECT_EQ(matrix.get_element(i, j), count++);
    }
  }
}

TEST(matrix_transformation, matrix_of_complements) {
  S21Matrix matrix = get_matrix_for_complements();

  EXPECT_EQ(get_matrix_of_compliments().EqMatrix(matrix.CalcComplements()),
            true);
}

TEST(matrix_transformation, matrix_of_complements_2) {
  S21Matrix matrix(4, 5);

  EXPECT_THROW(matrix.CalcComplements(), Array_exception);
}

TEST(matrix_transformation, matrix_inversion) {
  S21Matrix matrix = get_matrix_for_inversion();

  EXPECT_EQ(get_inverted_matrix().EqMatrix(matrix.InverseMatrix()), true);
}

TEST(matrix_transformation, matrix_inversion_2) {
  S21Matrix matrix(1, 1);
  matrix(0, 0) = 1.0;

  EXPECT_TRUE(matrix.EqMatrix(matrix.InverseMatrix()));
}

TEST(matrix_transformation, matrix_inversion_exceptions) {
  S21Matrix matrix(3, 3);
  matrix.get_element(0, 0) = 1.0;
  matrix.get_element(0, 1) = 2.0;
  matrix.get_element(0, 2) = 3.0;
  matrix.get_element(1, 0) = 4.0;
  matrix.get_element(1, 1) = 5.0;
  matrix.get_element(1, 2) = 6.0;
  matrix.get_element(2, 0) = 7.0;
  matrix.get_element(2, 1) = 8.0;
  matrix.get_element(2, 2) = 9.0;

  EXPECT_THROW(matrix.InverseMatrix(), Array_exception);
}

TEST(matrix_transformation, matrix_inversion_exceptions_2) {
  S21Matrix matrix(3, 4);

  EXPECT_THROW(matrix.InverseMatrix(), Array_exception);
}

TEST(matrix_operators, indexation_operator) {
  S21Matrix matrix(2, 2);
  matrix(0, 0) = 0.0;
  matrix(0, 1) = 1.0;
  matrix(1, 0) = 2.0;
  matrix(1, 1) = 3.0;

  EXPECT_EQ(matrix(0, 0), 0.0);
  EXPECT_EQ(matrix(0, 1), 1.0);
  EXPECT_EQ(matrix(1, 0), 2.0);
  EXPECT_EQ(matrix(1, 1), 3.0);

  const S21Matrix matrix2 = matrix;

  EXPECT_EQ(matrix2(0, 0), 0.0);
  EXPECT_EQ(matrix2(0, 1), 1.0);
  EXPECT_EQ(matrix2(1, 0), 2.0);
  EXPECT_EQ(matrix2(1, 1), 3.0);
}

TEST(matrix_operators, indexation_exceptions) {
  S21Matrix matrix(2, 2);
  matrix(0, 0) = 0.0;
  matrix(0, 1) = 1.0;
  matrix(1, 0) = 2.0;
  matrix(1, 1) = 3.0;

  EXPECT_THROW(matrix(-1, 1), Array_exception);
  EXPECT_THROW(matrix(3, 1), Array_exception);
  EXPECT_THROW(matrix(1, -1), Array_exception);
  EXPECT_THROW(matrix(1, 3), Array_exception);

  const S21Matrix matrix2 = matrix;

  EXPECT_THROW(matrix2(-1, 1), Array_exception);
  EXPECT_THROW(matrix2(3, 1), Array_exception);
  EXPECT_THROW(matrix2(1, -1), Array_exception);
  EXPECT_THROW(matrix2(1, 3), Array_exception);
}

TEST(matrix_operators, matrix_addition) {
  S21Matrix matrix1 = get_basic_matrix(3, 3);
  S21Matrix matrix2 = get_basic_matrix(3, 3);

  EXPECT_TRUE(get_basic_sum(3, 3).EqMatrix(matrix1 + matrix2));
}

TEST(matrix_operators, matrix_subtraction) {
  S21Matrix matrix1 = get_basic_matrix(3, 3);
  S21Matrix matrix2 = get_basic_matrix(3, 3);

  EXPECT_TRUE(get_basic_dif(3, 3).EqMatrix(matrix1 - matrix2));
}

TEST(matrix_operators, matrix_multiplication) {
  S21Matrix matrix1 = get_matrix_32();
  S21Matrix matrix2 = get_matrix_23();

  EXPECT_TRUE(get_basic_product().EqMatrix(matrix1 * matrix2));
}

TEST(matrix_operators, matrix_equality) {
  S21Matrix A(3, 3);
  S21Matrix B(3, 3);

  EXPECT_TRUE(A.EqMatrix(B));

  A.get_element(1, 1) = 1.0;
  EXPECT_FALSE(A == B);

  B.get_element(1, 1) = 1.000001;
  EXPECT_FALSE(A == B);

  B.get_element(1, 1) = 1.0000001;
  EXPECT_FALSE(A == B);

  B.get_element(1, 1) = 1.00000001;
  EXPECT_TRUE(A == B);
}

TEST(matrix_operators, addition_assignment) {
  S21Matrix matrix1 = get_basic_matrix(3, 3);
  S21Matrix matrix2 = get_basic_matrix(3, 3);

  matrix1 += matrix2;

  EXPECT_TRUE(get_basic_sum(3, 3).EqMatrix(matrix1));
}

TEST(matrix_operators, difference_assignment) {
  S21Matrix matrix1 = get_basic_matrix(3, 3);
  S21Matrix matrix2 = get_basic_matrix(3, 3);

  matrix1 -= matrix2;

  EXPECT_TRUE(get_basic_dif(3, 3).EqMatrix(matrix1));
}

TEST(matrix_operators, multiplication_matrix_assignment) {
  S21Matrix matrix1 = get_matrix_32();
  S21Matrix matrix2 = get_matrix_23();

  matrix1 *= matrix2;

  EXPECT_TRUE(get_basic_product().EqMatrix(matrix1));
}

TEST(matrix_operators, multiplication_number_assignment) {
  S21Matrix matrix = get_basic_matrix(3, 3);
  S21Matrix basic_matrix = get_basic_matrix(3, 3);

  matrix *= 3.14;

  for (int i = 0; i < matrix.get_rows(); ++i) {
    for (int j = 0; j < matrix.get_cols(); ++j) {
      EXPECT_EQ(matrix.get_element(i, j),
                basic_matrix.get_element(i, j) * 3.14);
    }
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}