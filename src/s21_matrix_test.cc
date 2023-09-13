#include <gtest/gtest.h>

#include <cstdlib>
#include <iostream>

#include "s21_matrix_oop.h"

double dum_sum_1_1[] = {-2,  -4,  -6,  -8,  -10, -12, -14, -16,
                        -18, -20, -22, -24, -0,  -0,  -0};

double dum_sum_1_2[] = {-1, -2,  -3,  -4,  -5, -6, -7, -8,
                        -9, -10, -11, -12, -0, -0, -0};

TEST(core_function_test, EqMatrix_1) {
  S21Matrix matrix_1, matrix_2;
  for (int i = 0; i < matrix_1.getRows(); i++) {
    for (int j = 0; j < matrix_1.getCols(); j++) {
      matrix_1(i, j) = 1;
      matrix_2(i, j) = 1;
    }
  }
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, EqMatrix_2) {
  S21Matrix matrix_1(5, 3), matrix_2(5, 3);
  matrix_1.InsertArray(dum_sum_1_1);
  matrix_2.InsertArray(dum_sum_1_2);
  EXPECT_FALSE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, EqMatrix_3) {
  S21Matrix matrix_1(5, 3), matrix_2(4, 3);
  matrix_1.InsertArray(dum_sum_1_1);
  for (int i = 0; i < matrix_2.getRows(); i++) {
    for (int j = 0; j < matrix_2.getCols(); j++) {
      matrix_2(i, j) = 1;
    }
  }
  EXPECT_FALSE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, SumMatrix_1) {
  S21Matrix matrix_1(3, 3), matrix_2(3, 3);
  for (int i = 0; i < matrix_2.getRows(); i++) {
    for (int j = 0; j < matrix_2.getCols(); j++) {
      matrix_2(i, j) = 1;
      matrix_1(i, j) = 2;
    }
  }
  matrix_2.SumMatrix(matrix_2);
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, SumMatrix_2) {
  S21Matrix matrix_1(5, 3), matrix_2(4, 3);
  EXPECT_THROW(matrix_1.SumMatrix(matrix_2), std::logic_error);
}

TEST(core_function_test, SubMatrix_1) {
  S21Matrix matrix_1(5, 3), matrix_2(4, 3);
  EXPECT_THROW(matrix_1.SubMatrix(matrix_2), std::logic_error);
}

TEST(core_function_test, SubMatrix_2) {
  S21Matrix matrix_1(3, 3), matrix_2(3, 3);
  for (int i = 0; i < matrix_2.getRows(); i++) {
    for (int j = 0; j < matrix_2.getCols(); j++) {
      matrix_2(i, j) = 1;
      matrix_1(i, j) = 2;
    }
  }
  matrix_1.SubMatrix(matrix_2);
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, SubMatrix_3) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  matrix_1(0, 0) = 2;
  matrix_1(0, 1) = 4;
  matrix_1(1, 0) = 6;
  matrix_1(1, 1) = 8;
  matrix_2.SubMatrix(matrix_1);
  matrix_1 *= -1;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, MulNumber_1) {
  S21Matrix matrix_1(3, 3), matrix_2(3, 3);
  for (int i = 0; i < matrix_2.getRows(); i++) {
    for (int j = 0; j < matrix_2.getCols(); j++) {
      matrix_1(i, j) = 3;
      matrix_2(i, j) = 6;
    }
  }
  matrix_1.MulNumber(2);
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, MulNumber_2) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  matrix_1(0, 0) = 2;
  matrix_1(0, 1) = 4;
  matrix_1(1, 0) = 6;
  matrix_1(1, 1) = 8;
  matrix_1.MulNumber(0);
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, MulMatrix_1) {
  S21Matrix marix_1(2, 3), mtrix_2(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      marix_1(i, j) = count;
      mtrix_2(j, i) = count;
      count++;
    }
  }
  marix_1.MulMatrix(mtrix_2);
  S21Matrix res(2, 2);
  res(0, 0) = 14;
  res(0, 1) = 32;
  res(1, 0) = 32;
  res(1, 1) = 77;
  EXPECT_TRUE(marix_1.EqMatrix(res));
}

TEST(core_function_test, MulMatrix_2) {
  S21Matrix matrix_1(2, 3), matrix_2(4, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      matrix_2(j, i) = count;
      count++;
    }
  }
  EXPECT_THROW(matrix_1.MulMatrix(matrix_2), std::logic_error);
}

TEST(core_function_test, Transpose_1) {
  S21Matrix matrix_1(4, 4), matrix_2(4, 4);
  double count = 1;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      matrix_1(i, j) = count;
      matrix_2(j, i) = count;
      count++;
    }
  }
  matrix_2 = matrix_2.Transpose();
  EXPECT_TRUE(matrix_2.EqMatrix(matrix_1));
}

TEST(core_function_test, Transpose_2) {
  S21Matrix matrix_1(2, 3), matrix_2(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      matrix_2(j, i) = count;
      count++;
    }
  }
  matrix_2 = matrix_2.Transpose();
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, CalcComplements_1) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  int count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  matrix_1 = matrix_1.CalcComplements();
  matrix_2(0, 0) = 4;
  matrix_2(0, 1) = -3;
  matrix_2(1, 0) = -2;
  matrix_2(1, 1) = 1;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, CalcComplements_2) {
  S21Matrix matrix_1(2, 3);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  EXPECT_THROW(matrix_1.CalcComplements(), std::logic_error);
}

TEST(core_function_test, Determinant_1) {
  S21Matrix matrix_1(3, 3);
  int count = 1;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  EXPECT_DOUBLE_EQ(0, matrix_1.Determinant());
}

TEST(core_function_test, Determinant_2) {
  S21Matrix matrix_1(3, 3);
  int count = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      count += i + j;
    }
    count *= 2;
  }
  EXPECT_DOUBLE_EQ(-12, matrix_1.Determinant());
}

TEST(core_function_test, Inverse_1) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  int count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  matrix_1 = matrix_1.InverseMatrix();
  matrix_2(0, 0) = -2;
  matrix_2(0, 1) = 1;
  matrix_2(1, 0) = 1.5;
  matrix_2(1, 1) = -0.5;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, Inverse_2) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  matrix_1(0, 0) = 1;
  matrix_1(0, 1) = 1;
  matrix_1(1, 0) = 3;
  matrix_1(1, 1) = 2;
  matrix_1 = matrix_1.InverseMatrix();
  matrix_2(0, 0) = -2;
  matrix_2(0, 1) = 1;
  matrix_2(1, 0) = 3;
  matrix_2(1, 1) = -1;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(get_set, t1) {
  S21Matrix m1(2, 2);
  m1(1, 1) = 123.456;
  EXPECT_EQ(m1.getCols(), 2);
  EXPECT_EQ(m1.getRows(), 2);
  EXPECT_DOUBLE_EQ(m1(1, 1), 123.456);
}

TEST(set_col, t1) {
  S21Matrix m(2, 2);
  m(1, 1) = 123.456;
  m.setCols(3);
  m(1, 2) = 999.999;
  EXPECT_EQ(m.getCols(), 3);
  EXPECT_DOUBLE_EQ(m(1, 1), 123.456);
  EXPECT_DOUBLE_EQ(m(1, 2), 999.999);
}

TEST(set_row, t1) {
  S21Matrix m1(2, 2);
  m1(1, 0) = 123.456;
  m1.setRows(3);
  m1(2, 0) = m1(1, 0);
  EXPECT_EQ(m1.getRows(), 3);
  EXPECT_EQ(m1.getCols(), 2);
  EXPECT_DOUBLE_EQ(m1(1, 0), m1(2, 0));
}

TEST(constructors, t1) {
  S21Matrix m1;
  EXPECT_EQ(m1.getCols(), 0);
  EXPECT_EQ(m1.getRows(), 0);

  S21Matrix m2(5, 5);
  m2(1, 1) = 123.456;
  EXPECT_EQ(m2.getCols(), 5);
  EXPECT_EQ(m2.getRows(), 5);

  S21Matrix m3(m2);
  EXPECT_EQ(m2 == m3, true);

  S21Matrix m4(std::move(m3));
  EXPECT_EQ(m3.getCols(), 0);
  EXPECT_EQ(m3.getRows(), 0);
  EXPECT_EQ(m4.getRows(), 5);
  EXPECT_EQ(m4.getCols(), 5);
  EXPECT_EQ(m4(1, 1), m2(1, 1));
}

TEST(func1, eq) {
  S21Matrix m1(123, 123);
  S21Matrix m2(m1);

  EXPECT_EQ(m1.EqMatrix(m2), true);
}

TEST(func2, summat) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);

  m1(0, 0) = 1;
  m1(0, 2) = 2;
  m1(2, 0) = 3;

  m2.SumMatrix(m1);

  EXPECT_EQ(m2(0, 0), 1);
  EXPECT_EQ(m2(0, 2), 2);
  EXPECT_EQ(m2(2, 0), 3);
}

TEST(func3, submat) {
  S21Matrix m1(3, 3);
  S21Matrix m2(3, 3);

  m1(0, 0) = 1;
  m1(0, 2) = 2;
  m1(2, 0) = 3;

  m2.SubMatrix(m1);

  EXPECT_EQ(m2(0, 0), -1);
  EXPECT_EQ(m2(0, 2), -2);
  EXPECT_EQ(m2(2, 0), -3);
}

TEST(func4, mulnum) {
  S21Matrix m2(3, 3);

  m2(0, 0) = 1.3;
  m2(0, 2) = 2.1;
  m2(2, 0) = 3.3;

  m2.MulNumber(0.412);

  EXPECT_EQ(m2(0, 0), 1.3 * 0.412);
  EXPECT_EQ(m2(0, 2), 2.1 * 0.412);
  EXPECT_EQ(m2(2, 0), 3.3 * 0.412);
}

TEST(func5, mulmat) {
  S21Matrix m1(2, 1);
  S21Matrix m2(1, 2);

  m1(0, 0) = 32;
  m1(1, 0) = 12;

  m2(0, 0) = 2;
  m2(0, 1) = 3;

  m1.MulMatrix(m2);

  EXPECT_DOUBLE_EQ(m1(0, 0), 32 * 2);
  EXPECT_DOUBLE_EQ(m1(1, 1), 12 * 3);
}

TEST(func6, transpose) {
  int rows = 2;
  int cols = 3;

  S21Matrix m(rows, cols);
  S21Matrix res(cols, rows);

  for (int i = 0, c = 1; i < rows; i++) {
    for (int j = 0; j < cols; j++, c++) {
      m(i, j) = c;
      res(j, i) = c;
    }
  }

  m = m.Transpose();

  ASSERT_TRUE(m == res);
}

TEST(func7, calccal) {
  S21Matrix m(3, 12);
  EXPECT_ANY_THROW(m.CalcComplements());
}

TEST(func7, calcalc) {
  int rows = 3;
  int cols = 3;

  S21Matrix given(rows, cols);
  S21Matrix expected(rows, cols);

  given(0, 0) = 1.0;
  given(0, 1) = 2.0;
  given(0, 2) = 3.0;
  given(1, 0) = 0.0;
  given(1, 1) = 4.0;
  given(1, 2) = 2.0;
  given(2, 0) = 5.0;
  given(2, 1) = 2.0;
  given(2, 2) = 1.0;

  expected(0, 0) = 0.0;
  expected(0, 1) = 10.0;
  expected(0, 2) = -20.0;
  expected(1, 0) = 4.0;
  expected(1, 1) = -14.0;
  expected(1, 2) = 8.0;
  expected(2, 0) = -8.0;
  expected(2, 1) = -2.0;
  expected(2, 2) = 4.0;

  S21Matrix res = given.CalcComplements();

  ASSERT_TRUE(res == expected);
}

TEST(func8, det) {
  S21Matrix m(5, 5);

  for (uint32_t i = 0; i < 5; i++)
    for (uint32_t j = 0; j < 5; j++) m(i, j) = j;

  ASSERT_NEAR(0, m.Determinant(), 1e-06);
}

TEST(func8, det1) {
  int size = 5;
  S21Matrix m(size, size);

  m(0, 1) = 6;
  m(0, 2) = -2;
  m(0, 3) = -1;
  m(0, 4) = 5;
  m(1, 3) = -9;
  m(1, 4) = -7;
  m(2, 1) = 15;
  m(2, 2) = 35;
  m(3, 1) = -1;
  m(3, 2) = -11;
  m(3, 3) = -2;
  m(3, 4) = 1;
  m(4, 0) = -2;
  m(4, 1) = -2;
  m(4, 2) = 3;
  m(4, 4) = -2;

  double res = m.Determinant();
  ASSERT_NEAR(res, 2480, 1e-6);
}

TEST(func8, det2) {
  int size = 3;
  S21Matrix m(size, size);

  m(0, 0) = 2;
  m(0, 1) = 3;
  m(0, 2) = 1;
  m(1, 0) = 7;
  m(1, 1) = 4;
  m(1, 2) = 1;
  m(2, 0) = 9;
  m(2, 1) = -2;
  m(2, 2) = 1;

  double res = m.Determinant();
  ASSERT_NEAR(res, -32, 1e-6);
}

TEST(func9, inverse) {
  int size = 3;
  S21Matrix given(size, size);
  S21Matrix expected(size, size);

  expected(0, 0) = 44300.0 / 367429.0;
  expected(0, 1) = -236300.0 / 367429.0;
  expected(0, 2) = 200360.0 / 367429.0;
  expected(1, 0) = 20600.0 / 367429.0;
  expected(1, 1) = 56000.0 / 367429.0;
  expected(1, 2) = -156483.0 / 367429.0;
  expected(2, 0) = 30900.0 / 367429.0;
  expected(2, 1) = 84000.0 / 367429.0;
  expected(2, 2) = -51010.0 / 367429.0;

  given(0, 0) = 2.8;
  given(0, 1) = 1.3;
  given(0, 2) = 7.01;
  given(1, 0) = -1.03;
  given(1, 1) = -2.3;
  given(1, 2) = 3.01;
  given(2, 0) = 0;
  given(2, 1) = -3;
  given(2, 2) = 2;

  ASSERT_TRUE(given.InverseMatrix() == expected);
}

TEST(operator_overloading, operator_sum) {
  S21Matrix matrix_1(3, 3), matrix_2(3, 3);
  for (int i = 0; i < matrix_2.getRows(); i++) {
    for (int j = 0; j < matrix_2.getCols(); j++) {
      matrix_2(i, j) = 1;
      matrix_1(i, j) = 2;
    }
  }
  matrix_2 += matrix_2;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(operator_overloading, operator_minus) {
  S21Matrix matrix_1(3, 3), matrix_2(3, 3);
  for (int i = 0; i < matrix_2.getRows(); i++) {
    for (int j = 0; j < matrix_2.getCols(); j++) {
      matrix_2(i, j) = 1;
      matrix_1(i, j) = 2;
    }
  }
  matrix_1 -= matrix_2;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(operator_overloading, operator_mult) {
  S21Matrix matrix_1(2, 3), matrix_2(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      matrix_2(j, i) = count;
      count++;
    }
  }
  matrix_1 *= matrix_2;
  S21Matrix res(2, 2);
  res(0, 0) = 14;
  res(0, 1) = 32;
  res(1, 0) = 32;
  res(1, 1) = 77;
  EXPECT_TRUE(matrix_1.EqMatrix(res));
}

TEST(operator_overloading, operator_copy) {
  S21Matrix matrix_1(2, 3), matrix_2(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  matrix_2 = matrix_1;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(operator_overloading, operator_sum_both) {
  S21Matrix matrix_1(3, 3), matrix_2(3, 3);
  for (int i = 0; i < matrix_2.getRows(); i++) {
    for (int j = 0; j < matrix_2.getCols(); j++) {
      matrix_2(i, j) = 1;
      matrix_1(i, j) = 2;
    }
  }
  matrix_2 = matrix_2 + matrix_2;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(operator_overloading, operator_sub_both) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  matrix_1(0, 0) = 2;
  matrix_1(0, 1) = 4;
  matrix_1(1, 0) = 6;
  matrix_1(1, 1) = 8;
  matrix_2 = matrix_2 - matrix_1;
  matrix_1 *= -1;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(operator_overloading, operator_mult_both) {
  S21Matrix matrix_1(2, 3), matrix_2(3, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      matrix_2(j, i) = count;
      count++;
    }
  }
  matrix_1 = matrix_1 * matrix_2;
  S21Matrix res(2, 2);
  res(0, 0) = 14;
  res(0, 1) = 32;
  res(1, 0) = 32;
  res(1, 1) = 77;
  EXPECT_TRUE(matrix_1.EqMatrix(res));
}

TEST(operator_overloading, operator_mult_number) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  matrix_1 = matrix_1 * 2;
  matrix_2(0, 0) = 2;
  matrix_2(0, 1) = 4;
  matrix_2(1, 0) = 6;
  matrix_2(1, 1) = 8;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(throw_test, S21Matrix_throw) {
  EXPECT_THROW(S21Matrix matrix_1(-2, -1), std::logic_error);
}

TEST(throw_test, Determinant_throw) {
  S21Matrix matrix_1(2, 3);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  EXPECT_THROW(matrix_1.Determinant(), std::logic_error);
}

TEST(throw_test, InverseMatrix_throw) {
  S21Matrix matrix_1(3, 3);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      if (i == 0) {
        matrix_1(i, j) = 0;
      } else
        matrix_1(i, j) = count;
      count++;
    }
  }
  EXPECT_THROW(matrix_1.InverseMatrix(), std::logic_error);
}

TEST(throw_test, InverseMatrix_throw_1) {
  S21Matrix matrix_1(1, 1), matrix_2(1, 1);
  double count = 1;
  for (int i = 0; i < 1; i++) {
    for (int j = 0; j < 1; j++) {
      matrix_1(i, j) = 1;
      matrix_2(i, j) = 2;
      count++;
    }
  }
  matrix_2 = matrix_1.InverseMatrix();
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(core_function_test, setRows_1) {
  S21Matrix matrix_1(3, 3);
  EXPECT_THROW(matrix_1.setRows(-1), std::logic_error);
}

TEST(core_function_test, setCols) {
  S21Matrix matrix_1(3, 3);
  EXPECT_THROW(matrix_1.setCols(-1), std::logic_error);
}

TEST(core_function_test, operator_1) {
  S21Matrix matrix_1(3, 3);
  EXPECT_THROW(matrix_1(-1, -2), std::logic_error);
}

TEST(core_function_test, operator_2) {
  S21Matrix matrix_1(3, 3);
  EXPECT_THROW(matrix_1(4, 5), std::logic_error);
}

TEST(addition, big_number_test) {
  S21Matrix matrix_check(999, 1999);
  ASSERT_NEAR(matrix_check(998, 1998), 0, 1e-7);
}

TEST(throw_test, mul_number_matrix_test) {
  S21Matrix matrix_1(2, 2), matrix_2(2, 2);
  double count = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      matrix_1(i, j) = count;
      count++;
    }
  }
  matrix_1 = 2 * matrix_1;
  matrix_2(0, 0) = 2;
  matrix_2(0, 1) = 4;
  matrix_2(1, 0) = 6;
  matrix_2(1, 1) = 8;
  EXPECT_TRUE(matrix_1.EqMatrix(matrix_2));
}

TEST(operator_overloading, operator_index) {
  const S21Matrix matrix_1(2, 2);
  const double x = matrix_1(1, 0);
  EXPECT_TRUE(x == 0);
  EXPECT_THROW(matrix_1(3, 0), std::logic_error);
  EXPECT_ANY_THROW(matrix_1(0, -3));
}

TEST(sega, operator_index) {
  S21Matrix matrix_1;
  EXPECT_THROW(matrix_1(0, 0), std::logic_error);
}

int main() {
  testing::InitGoogleTest();
  int retval = RUN_ALL_TESTS();
  return retval;
}
