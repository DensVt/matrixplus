#include "s21_matrix_oop.h"

// Конструкторы
S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}
S21Matrix::S21Matrix(long int rows, long int cols) : rows_(rows), cols_(cols) {
  // Проверка на корректность размеров
  if (rows_ <= 0 || cols_ <= 0) {
    throw std::out_of_range("Matrix has no value or size less than 1X1");
  }
  // Вызов метода для создания матрицы
  CreateMatrix();
}

// Копирование и перемещение
S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  CopyMatrix(other);
}
S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

// Деструктор
S21Matrix::~S21Matrix() {
  if (matrix_) {
    RemoveMatrix(*this);
  }
}

// Создание матрицы в динамической памяти
void S21Matrix::CreateMatrix() {
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    try {
      matrix_[i] = new double[cols_]();
    } catch (...) {
      for (int j = 0; j < i; j++) delete[] matrix_[j];
      delete[] matrix_;
      throw;
    }
  }
}

// Удаление матрицы, освобождение памяти
void S21Matrix::RemoveMatrix(S21Matrix& other) noexcept {
  for (int i = 0; i < other.rows_; i++) {
    delete[] other.matrix_[i];
  }
  delete[] other.matrix_;
}

// Копирование матрицы
void S21Matrix::CopyMatrix(const S21Matrix& other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  CreateMatrix();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

// Сравнение матриц на равенство
bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  if (!this->EqualColRow(other)) {
    return false;
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (std::abs(matrix_[i][j] - other.matrix_[i][j]) >= min_diff_) {
        return false;
      }
    }
  }
  return true;
}

// Сложение матриц
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (!this->EqualColRow(other)) {
    throw std::out_of_range("Error! Matrixes should have the same size");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

// Вычитание матриц
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (!EqualColRow(other)) {
    throw std::out_of_range("Error! Matrixes should have the same size");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
    }
  }
}

// Умножение матрицы на число
void S21Matrix::MulNumber(const double num) const noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

// Проверка равенства размеров матриц
bool S21Matrix::EqualColRow(const S21Matrix& other) const noexcept {
  return other.cols_ == cols_ && other.rows_ == rows_;
}

// Умножение матриц
void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->cols_ != other.rows_) {
    throw std::out_of_range(
        "Matrixes should have same size (matrix_1:col == matrix_2:rows_");
  }
  S21Matrix buff(this->rows_, other.cols_);
  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int q = 0; q < other.rows_; q++) {
        buff.matrix_[i][j] += matrix_[i][q] * other.matrix_[q][j];
      }
    }
  }
  RemoveMatrix(*this);
  this->CopyMatrix(buff);
}

// Транспонирование матрицы
S21Matrix S21Matrix::Transpose() const {
  S21Matrix new_matrix(cols_, rows_);
  for (int i = 0; i < cols_; i++) {
    for (int j = 0; j < rows_; j++) {
      new_matrix.matrix_[i][j] = matrix_[j][i];
    }
  }
  return new_matrix;
}

// Расчет алгебраических дополнений для матрицы
S21Matrix S21Matrix::CalcComplements() const {
  if (this->rows_ != this->cols_) {
    throw std::out_of_range("Error! The matrix is not quadratic.");
  }
  S21Matrix result(rows_, cols_);
  int sign = 1;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      double buff;
      S21Matrix temp = GetCutMatrix(i, j);
      buff = temp.Determinant();
      result.matrix_[i][j] = buff * sign;
      sign = sign * -1;
    }
    if (this->rows_ % 2 == 0) {
      sign = sign * -1;
    }
  }
  return result;
}

// Вычисление определителя матрицы
double S21Matrix::Determinant() const {
  if (cols_ != rows_) {
    throw std::out_of_range("Error! The matrix is not quadratic.");
  } else if (cols_ == 1) {
    return matrix_[0][0];
  } else {
    double result = 0;
    int sign = 1;
    for (int j = 0; j < cols_; j++) {
      double buff = 0;
      S21Matrix temp = GetCutMatrix(0, j);
      buff = temp.Determinant();
      buff = buff * matrix_[0][j] * sign;
      result += buff;
      sign = sign * -1;
    }
    return result;
  }
}

// Вычисление обратной матрицы
S21Matrix S21Matrix::InverseMatrix() const {
  double determ = this->Determinant();
  if (std::abs(determ) < min_diff_) {
    throw std::out_of_range("Determination = 0");
  }
  S21Matrix buff;
  if (rows_ == 1 || cols_ == 1) {
    buff.CopyMatrix(*this);
  } else {
    buff = Transpose();
    buff = buff.CalcComplements();
  }
  buff.MulNumber(1 / determ);
  return buff;
}

// Установка нового числа строк
void S21Matrix::setRows(int len_row) {
  if (len_row <= 0) {
    throw std::out_of_range("Error, size of len_row should be positive");
  }
  if (len_row != rows_) {
    double** buff_matrix = new double*[len_row];
    for (int i = 0; i < len_row; i++) {
      buff_matrix[i] = new double[cols_];
      for (int j = 0; j < cols_; j++) {
        if (i < rows_) {
          buff_matrix[i][j] = matrix_[i][j];
        } else {
          buff_matrix[i][j] = 0;
        }
      }
    }
    RemoveMatrix(*this);
    rows_ = len_row;
    matrix_ = buff_matrix;
  }
}

// Установка нового числа столбцов
void S21Matrix::setCols(int len_col) {
  if (len_col <= 0) {
    throw std::out_of_range("Error, size of len_col should be positive");
  }
  if (len_col != cols_) {
    double** buff_matrix = new double*[rows_];
    for (int i = 0; i < this->rows_; i++) {
      buff_matrix[i] = new double[len_col];
      for (int j = 0; j < len_col; j++) {
        if (j < cols_) {
          buff_matrix[i][j] = matrix_[i][j];
        } else {
          buff_matrix[i][j] = 0;
        }
      }
    }
    RemoveMatrix(*this);
    cols_ = len_col;
    matrix_ = buff_matrix;
  }
}

int S21Matrix::getRows() const noexcept {
  return this->rows_;
}  // Получение числа строк
int S21Matrix::getCols() const noexcept {
  return this->cols_;
}  // Получение числа столбцов

// Перегрузка операторов для работы с матрицами
S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double number) const noexcept {
  S21Matrix result(*this);
  result.MulNumber(number);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix operator*(const double number, const S21Matrix& other) noexcept {
  S21Matrix result(other);
  result.PrintMatrix();
  result.MulNumber(number);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) const noexcept {
  return this->EqMatrix(other);
}

// Перегрузка оператора присваивания
S21Matrix& S21Matrix::operator=(S21Matrix& other) {
  if (this != &other) {
    RemoveMatrix(*this);
    CopyMatrix(other);
  }
  return *this;
}

// Перегрузка оператора присваивания для перемещения
S21Matrix& S21Matrix::operator=(S21Matrix&& other) noexcept {
  if (this != &other) {
    RemoveMatrix(*this);
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;

    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }
  return *this;
}

// Перегрузка оператора += для суммирования матриц
S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

// Перегрузка оператора -= для вычитания матриц
S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

// Перегрузка оператора *= для умножения матриц
S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

// Перегрузка оператора *= для умножения матрицы на число
S21Matrix& S21Matrix::operator*=(const double number) noexcept {
  this->MulNumber(number);
  return *this;
}

// Оператор доступа к элементам матрицы (для изменения)
double& S21Matrix::operator()(const int i, const int j) {
  if (rows_ <= i || cols_ <= j) {
    throw std::out_of_range("Error! Index is out of range.");
  } else if (i < 0 || j < 0) {
    throw std::out_of_range("Error! Index is lower than 0.");
  }
  return matrix_[i][j];
}

// Оператор доступа к элементам матрицы (для чтения)
const double& S21Matrix::operator()(const int i, const int j) const {
  if (rows_ <= i || cols_ <= j) {
    throw std::out_of_range("Error! Index is out of range.");
  } else if (i < 0 || j < 0) {
    throw std::out_of_range("Error! Index is lower than 0.");
  }
  return matrix_[i][j];
}

// Вставка массива в матрицу
void S21Matrix::InsertArray(double* dum) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = dum[i * cols_ + j];
    }
  }
}

// Вывод матрицы на экран
void S21Matrix::PrintMatrix() {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      std::cout << matrix_[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

// Получение матрицы, исключая определенную строку и столбец
S21Matrix S21Matrix::GetCutMatrix(int row, int col) const {
  int buff_row = 0;
  S21Matrix buff_matrix(rows_ - 1, cols_ - 1);
  for (int i = 0; i < rows_; i++) {
    int column = 0;
    if (i != row) {
      for (int j = 0; j < cols_; j++) {
        if (j != col) {
          buff_matrix.matrix_[buff_row][column++] = matrix_[i][j];
        }
      }
      buff_row++;
    }
  }
  return buff_matrix;
}