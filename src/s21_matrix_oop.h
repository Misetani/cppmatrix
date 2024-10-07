#pragma once

class S21Matrix {
 private:
  int rows = 0;
  int cols = 0;
  double** matrix = nullptr;

 public:
  /* ======================================================================= */
  /*                      Constructors and Destructors                       */
  /* ======================================================================= */

  /**
   * @brief A basic constructor. Default size is 0x0
   *
   * @note No need to define it explicitly, but this is part of the task.
   */
  S21Matrix() {}

  /**
   * @brief Parametrized constructor with number of rows and columns
   */
  S21Matrix(int r, int c);

  /**
   * @brief Copy constructor
   */
  S21Matrix(const S21Matrix& other);

  /**
   * @brief Move constructor
   */
  S21Matrix(S21Matrix&& other);

  /**
   * @brief Destructor
   */
  ~S21Matrix();

  /* ======================================================================= */
  /*                           Getters and Setters                           */
  /* ======================================================================= */

  int get_rows() const { return rows; }

  int get_cols() const { return cols; }

  /**
   * @brief Sets the number of rows
   *
   * @note If the matrix increases in size, it is filled with zeros.
   * If it decreases in size, the excess is simply discarded;
   */
  void set_rows(int r);

  /**
   * @brief Sets the number of columns
   *
   * @note If the matrix increases in size, it is filled with zeros.
   * If it decreases in size, the excess is simply discarded;
   */
  void set_cols(int c);

  /* ======================================================================= */
  /*                              Main Methods                               */
  /* ======================================================================= */

  /**
   * @brief Checks matrices for equality with each other.
   */
  bool EqMatrix(const S21Matrix& other) const;

  /**
   * @brief Adds the second matrix to the current one
   * @throws Array_exception: different matrix dimensions
   */
  void SumMatrix(const S21Matrix& other);

  /**
   * @brief Subtracts another matrix from the current one
   * @throws Array_exception: different matrix dimensions
   */
  void SubMatrix(const S21Matrix& other);

  /**
   * @brief Multiplies the current matrix by a number
   */
  void MulNumber(const double num);

  /**
   * @brief Multiplies the current matrix by another matrix
   * @throws Array_exception: The number of columns of the first matrix is
   * not equal to the number of rows of the second matrix
   */
  void MulMatrix(const S21Matrix& other);

  /**
   * @brief Creates a new transposed matrix from the current one and returns it
   */
  S21Matrix Transpose() const;

  /**
   * @brief Calculates the algebraic addition matrix of the current one and
   * returns it
   * @throws Array_exception: The matrix is not square
   */
  S21Matrix CalcComplements() const;

  /**
   * @brief Calculates and returns the determinant of the current matrix
   * @throws Array_exception: The matrix is not square
   */
  double Determinant() const;

  /**
   * @brief Calculates and returns the inverse matrix
   * @throws Array_exception: Matrix determinant is 0
   */
  S21Matrix InverseMatrix() const;

  /* ======================================================================= */
  /*                            Matrix Operators                             */
  /* ======================================================================= */

  /**
   * @brief Addition of two matrices.
   * @throws Array_exception: Different matrix dimensions.
   */
  S21Matrix operator+(const S21Matrix& other) const;

  /**
   * @brief Subtraction of one matrix from another.
   * @throws Array_exeption: Different matrix dimensions.
   */
  S21Matrix operator-(const S21Matrix& other) const;

  /**
   * @brief Matrix multiplication and matrix multiplication by a number.
   * @throws Array_exception: The number of columns of the first matrix does
   * not equal the number of rows of the second matrix.
   */
  S21Matrix operator*(const S21Matrix& other) const;

  /**
   * @brief Checks for matrices equality (EqMatrix).
   */
  bool operator==(const S21Matrix& other) const;

  /**
   * @brief Performs deep copying from one matrix to another
   *
   * @param other Right operand
   * @returns Left operand
   */
  S21Matrix& operator=(const S21Matrix& other);

  /**
   * @brief Addition assignment (SumMatrix)
   * @throws Array_exception: different matrix dimensions.
   */
  S21Matrix& operator+=(const S21Matrix& other);

  /**
   * @brief Difference assignment (SubMatrix)
   * @throws Array_exception: different matrix dimensions.
   */
  S21Matrix& operator-=(const S21Matrix& other);

  /**
   * @brief Multiplication assignment (MulMatrix).
   * @throws Array_exception: The number of columns of the first matrix does
   * not equal the number of rows of the second matrix.
   */
  S21Matrix& operator*=(const S21Matrix& other);

  /**
   * @brief Multiplication assignment (MulNumber).
   * @throws Array_exception: The number of columns of the first matrix does
   * not equal the number of rows of the second matrix.
   */
  S21Matrix& operator*=(double num);

  /**
   * @brief Indexation by matrix elements (row, column).
   *
   * @throws Array_exception: Index is outside the matrix.
   */
  double& operator()(int i, int j) { return get_element(i, j); }

  /**
   * @brief Indexation by matrix elements (row, column).
   *
   * @throws Array_exception: Index is outside the matrix.
   * @note Const version
   */
  const double& operator()(int i, int j) const { return get_element(i, j); }

  /* ======================================================================= */
  /*                         Helper public methods                           */
  /* ======================================================================= */

  void print() const;

  /**
   * @brief Indexation by matrix elements (row, column).
   *
   * @throws Array_exception: Index is outside the matrix.
   */
  double& get_element(int row, int col);

  /**
   * @brief Indexation by matrix elements (row, column).
   *
   * @throws Array_exception: Index is outside the matrix.
   * @note Const version
   */
  const double& get_element(int row, int col) const;

 private:
  /* ======================================================================= */
  /*                         Helper private methods                          */
  /* ======================================================================= */

  void free();

  void init_matrix(int rows, int cols);
};