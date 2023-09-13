#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H
#include <iostream>

class S21Matrix {
 private:
  // Приватные переменные
  int rows_, cols_;  // Количество строк и столбцов
  double** matrix_;  // Указатель на двумерную матрицу
  const float min_diff_ =
      0.0000001;  // Минимальное значение для сравнения вещественных чисел

  // Приватные вспомогательные методы
  void CreateMatrix();  // Инициализация матрицы
  void CopyMatrix(const S21Matrix& other);  // Копирование матрицы
  void RemoveMatrix(S21Matrix& other) noexcept;  // Удаление матрицы

 public:
  // Конструкторы и деструктор
  S21Matrix() noexcept;
  S21Matrix(long int rows, long int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  // Accessors и Mutators
  int getRows() const noexcept;
  int getCols() const noexcept;
  void setRows(int rows);
  void setCols(int cols);
  bool EqualColRow(const S21Matrix& other) const noexcept;

  // Пользовательские методы
  bool EqMatrix(const S21Matrix& other) const noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num) const noexcept;
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;

  // Перегрузка операторов
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const double number) const noexcept;
  S21Matrix operator*(const S21Matrix& other) const;
  bool operator==(const S21Matrix& other) const noexcept;
  S21Matrix& operator=(S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double number) noexcept;
  double& operator()(const int i, const int j);
  const double& operator()(const int i, const int j) const;

  // Вспомогательные методы
  void InsertArray(double* dum);  // Вставка массива в матрицу
  void PrintMatrix();             // Вывод матрицы на экран
  S21Matrix GetCutMatrix(int row,
                         int col) const;  // Получение "вырезанной" матрицы
};

// Внешние операторы
S21Matrix operator*(
    const double number,
    const S21Matrix& other) noexcept;  // Умножение числа на матрицу
#endif                                 // S21_MATRIX_OOP_H