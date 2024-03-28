//
// Created by Covald on 19.03.2024.
//

#include "Matrix.hpp"
#include <stdexcept>
#include <cmath>
#include <iomanip>

double &Matrix::at(size_t row, size_t col) {
    return n_data[row * n_col + col];
}

double &Matrix::at(size_t row, size_t col) const {
    return n_data[row * n_col + col];
}

bool Matrix::init(size_t rows, size_t cols) {
    if (rows != 0 && cols != 0) {
        n_row = rows;
        n_col = cols;
        n_data = new double[rows * cols];

        setZero();
    }
    return is_valid;
}

Matrix::Matrix(size_t cols) {
    init(1, cols);
}

Matrix::Matrix(size_t rows, size_t cols) {
    init(rows, cols);
}

Matrix::~Matrix() {
    delete[] n_data;
}

Matrix::Matrix(const Matrix &mat) {
    if (!mat.isValid()) Matrix();
    bool is_created = init(mat.rows(), mat.cols());
    if (is_created)
        for (size_t row = 0; row < n_row; row++)
            for (size_t col = 0; col < n_col; col++)
                at(row, col) = mat.coeffRef(row, col);
}

Matrix Matrix::operator*(const Matrix &mat) const {
    if (n_row != mat.cols() || !is_valid || !mat.is_valid) return {};

    Matrix result(n_row, mat.cols());

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < mat.cols(); col++)
            for (size_t k = 0; k < n_col; k++)
                result.coeffRef(row, col) += at(row, k) * mat.coeffRef(k, col);

    return result;
}

Matrix Matrix::operator-(const Matrix &mat) const {
    if (n_row != mat.n_row || n_col != mat.n_col || !is_valid || !mat.is_valid) return {};

    Matrix result = *this;
    if (!result.data()) return {};

    result -= mat;
    return result;
}

Matrix Matrix::operator+(const Matrix &mat) const {
    if (n_row != mat.n_row || n_col != mat.n_col || !is_valid || !mat.is_valid) return {};

    Matrix result = *this;
    if (!result.data()) return {};

    result += mat;
    return result;
}

Matrix Matrix::operator*(double value) const {
    if (!is_valid) return {};

    Matrix result = *this;
    result *= value;
    return result;
}

Matrix Matrix::operator/(double value) const {
    if (!is_valid) return {};

    Matrix result = *this;
    result /= value;
    return result;
}

Matrix &Matrix::operator=(const Matrix &mat) {
    if (this == &mat) return *this;
    if (!mat.isValid()) throw std::runtime_error("Try to move non valid matrix.");
    delete[] n_data;
    n_row = mat.n_row;
    n_col = mat.n_col;
    n_data = new double[n_row * n_col];

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < n_col; col++)
            at(row, col) = mat.coeffRef(row, col);

    this->is_valid = true;
    return *this;
}

Matrix &Matrix::operator*=(const Matrix &mat) {
    // Решил таки делать нормально, ибо все равно несколько раз копировать, чтобы эту же матрицу ресайзнуть
    *this = *this * mat;
    return *this;
}

Matrix &Matrix::operator+=(const Matrix &mat) {
    if (!is_valid || !mat.isValid() || n_row != mat.rows() || n_col != mat.cols()) {
        is_valid = false;
        return *this;
    };

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < n_col; col++)
            at(row, col) += mat.coeffRef(row, col);

    return *this;
}

Matrix &Matrix::operator-=(const Matrix &mat) {
    if (!is_valid || !mat.isValid() || n_row != mat.rows() || n_col != mat.cols()) {
        is_valid = false;
        return *this;
    };

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < n_col; col++)
            at(row, col) -= mat.coeffRef(row, col);

    return *this;
}

Matrix &Matrix::operator*=(double value) {
    if (!is_valid) return *this;

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < n_col; col++)
            at(row, col) *= value;

    return *this;
}

Matrix &Matrix::operator/=(double value) {
    if (!is_valid) {
        return *this;
    };

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < n_col; col++)
            at(row, col) /= value;

    return *this;
}

bool Matrix::isValid() const {
    return is_valid;
}

void Matrix::resize(size_t rows, size_t cols) {
    delete[] n_data;
    init(rows, cols);
}

const double &Matrix::coeffRef(size_t rowIdx, size_t colIdx) const {
    if (rowIdx > n_row - 1) throw std::out_of_range("Row index out of range");
    if (colIdx > n_col - 1) throw std::out_of_range("Col index out of range");

    return n_data[rowIdx * n_col + colIdx];
}

double &Matrix::coeffRef(size_t rowIdx, size_t colIdx) {
    if (rowIdx > n_row - 1) throw std::out_of_range("Row index out of range");
    if (colIdx > n_col - 1) throw std::out_of_range("Col index out of range");

    return n_data[rowIdx * n_col + colIdx];
}

const double *Matrix::data() const {
    return n_data;
}

double *Matrix::data() {
    return n_data;
}

size_t Matrix::rows() const {
    return n_row;
}

size_t Matrix::cols() const {
    return n_col;
}

Matrix &Matrix::setIdentity() {
    if (n_row != n_col) return *this;
    if (!n_data) return *this;


    for (size_t row = 0; row < n_row; row++) {
        for (size_t col = 0; col < n_col; col++) {
            if (row == col) {
                n_data[row * n_col + col] = 1;
            }
            std::cout << n_data[row * n_col + col];
        }
        std::cout << std::endl;
    }

    return *this;
}

Matrix &Matrix::setZero() {
    if (!n_data) return *this;
//    std::fill_n(n_data, n_row * n_col, 0);

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < n_col; col++)
            at(row, col) = 0;

//    std::cout << "SetZero" << std::endl << *this;
    is_valid = true;
    return *this;
}

Matrix &Matrix::setConstants(double value) {
    if (!n_data) return *this;

    for (size_t row = 0; row < n_row; row++) {
        for (size_t col = 0; col < n_col; col++) {
            at(row, col) = value;
        }
    }

    return *this;
}

Matrix &Matrix::setIdentity(size_t rows, size_t cols) {
    if (rows != cols) throw std::invalid_argument("Rowe must be equal cols.");
    this->resize(rows, cols);
    this->setIdentity();
    return *this;
}

Matrix &Matrix::setZero(size_t rows, size_t cols) {
    this->resize(rows, cols);
    this->setZero();
    return *this;
}

Matrix &Matrix::setConstants(size_t rows, size_t cols, double value) {
    this->resize(rows, cols);
    this->setConstants(value);
    return *this;
}

Matrix Matrix::transpose() const {
    if (!is_valid) return {};
    Matrix temp(n_row, n_col);

    // if (!temp.isValid()) throw std::runtime_error("Cannot allocate memory for new matrix");

    for (size_t row = 0; row < n_row; row++) {
        for (size_t col = 0; col < n_col; col++) {
            temp.at(row, col) = coeffRef(col, row);
        }
    }
    return temp;
}

static double determinant(double **Arr, unsigned size) {
    unsigned i, j;
    double det = 0;
    double **matr;
    if (size == 1) {
        det = Arr[0][0];
    } else if (size == 2) {
        det = Arr[0][0] * Arr[1][1] - Arr[0][1] * Arr[1][0];
    } else {
        matr = new double *[size - 1];
        for (i = 0; i < size; ++i) {
            for (j = 0; j < size - 1; ++j) {
                if (j < i)
                    matr[j] = Arr[j];
                else
                    matr[j] = Arr[j + 1];
            }
            det += pow((double) -1, (i + j)) * determinant(matr, size - 1) * Arr[i][size - 1];
        }
        delete[] matr;
    }
    return det;
}

static void inversion(double **A, int N) {
    double temp;

    auto **E = new double *[N];

    for (int i = 0; i < N; i++)
        E[i] = new double[N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++) {
        temp = A[k][k];

        for (int j = 0; j < N; j++) {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++) {
            temp = A[i][k];

            for (int j = 0; j < N; j++) {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--) {
        for (int i = k - 1; i >= 0; i--) {
            temp = A[i][k];

            for (int j = 0; j < N; j++) {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete[] E[i];

    delete[] E;
}

Matrix Matrix::inverse() const {
    if (!is_valid || n_col != n_row) return {};

    size_t size = n_row;
    auto **temp = new double *[size];
    for (size_t row = 0; row < size; row++) {
        temp[row] = new double[size];
        for (size_t col = 0; col < size; col++) {
            temp[row][col] = at(row, col);
        }
    }

    inversion(temp, size);

    Matrix result = *this;

    for (size_t row = 0; row < n_row; row++)
        for (size_t col = 0; col < n_col; col++)
            result.at(row, col) = temp[row][col];

    for (size_t row = 0; row < size; row++)
        delete [] temp[row];

    delete [] temp;

    return result;
}

double determinant(double **matrix, size_t size) { // ToDo
    double det = 0;
    int sign = 1;

    // Base Case
    if (size == 1) {
        det = matrix[0][0];
    } else if (size == 2) {
        det = (matrix[0][0] * matrix[1][1])
              - (matrix[0][1] * matrix[1][0]);
    }

        // Perform the Laplace Expansion
    else {
        for (size_t i = 0; i < size; i++) {

            // Stores the cofactor matrix
            auto **cofactor = new double *[size - 1];
            for (size_t j = 0; j < size - 1; j++) {
                cofactor[j] = new double[size - 1];
            }
            int sub_i = 0, sub_j = 0;
            for (size_t j = 1; j < size; j++) {
                for (size_t k = 0; k < size; k++) {
                    if (k == i) {
                        continue;
                    }
                    cofactor[sub_i][sub_j] = matrix[j][k];
                    sub_j++;
                }
                sub_i++;
                sub_j = 0;
            }

            // Update the determinant value
            det += sign * matrix[0][i]
                   * determinant(cofactor, size - 1);
            sign = -sign;
            for (size_t j = 0; j < size - 1; j++) {
                delete[] cofactor[j];
            }
            delete[] cofactor;
        }
    }

    // Return the final determinant value
    return det;
}

double Matrix::det() const {
    if (!is_valid || n_col != n_row) return std::nan("NaN");
    size_t size = n_col;

    auto **temp = new double *[size];
    for (size_t row = 0; row < size; row++) {
        temp[row] = new double[size];
        for (size_t col = 0; col < size; col++) {
            temp[row][col] = at(row, col);
        }
    }

    double _det = determinant(temp, size);

    for (size_t row = 0; row < size; row++) {
        if (temp[row]) delete[] temp[row];
    }
    delete[] temp;

    return _det;
}

Matrix Matrix::identity(size_t rows, size_t cols) {
    if (rows != cols) throw std::invalid_argument("Rowe must be equal cols.");
    Matrix result = Matrix().setIdentity(rows, cols);
    return result;
}

Matrix Matrix::zeros(size_t rows, size_t cols) {
    return Matrix().setZero(rows, cols);
}

Matrix Matrix::constants(size_t rows, size_t cols, double value) {
    return Matrix(rows, cols).setConstants(value);
}

Matrix operator*(double value, const Matrix &mat) {
    return mat * value;
}

std::ostream &operator<<(std::ostream &stream, const Matrix &matrix) {
    for (size_t row = 0; row < matrix.rows(); row++) {
        stream << "|";
        for (size_t col = 0; col < matrix.cols(); col++) {
            stream << std::setw(2) << matrix.coeffRef(row, col) << std::setw(2);
            stream << "|";
        }
        stream << std::endl;
    }
    return stream;
}




