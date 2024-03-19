//
// Created by Covald on 19.03.2024.
//

#include "Matrix.hpp"
#include <stdexcept>
#include <cmath>
#include <cstring>


Matrix::Matrix() : nrow(0), ncol(0), is_valid(false), ndata(nullptr) {}

Matrix::Matrix(size_t cols) {
    if (cols == 0) {
        this->nrow = 0;
        this->ncol = 0;
        this->ndata = nullptr;
    } else {
        this->nrow = 1;
        this->ncol = cols;
        this->ndata = new double[this->nrow * this->ncol];
    }
    this->is_valid = false;
}

Matrix::Matrix(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) {
        this->nrow = 0;
        this->ncol = 0;
        this->ndata = nullptr;
    } else {
        this->nrow = rows;
        this->ncol = cols;
        this->ndata = new double[this->nrow * this->ncol];
    }
    this->is_valid = false;
}

Matrix::~Matrix() {
    delete[] this->ndata;
}

Matrix::Matrix(const Matrix &mat) {
    this->nrow = mat.nrow;
    this->ncol = mat.ncol;
    this->ndata = new double[this->nrow * this->ncol];

    std::memcpy(this->ndata, mat.ndata, this->nrow * this->ncol);

    this->is_valid = true;
}

Matrix Matrix::operator*(const Matrix &mat) {
    if (this->nrow != mat.ncol || !this->is_valid || !mat.is_valid) return {};

    auto *result = new Matrix(this->nrow, mat.ncol);
    if (!result->ndata) return {};

    for (size_t row = 0; row < result->nrow; row++) {
        for (size_t col = 0; col < result->ncol; col++) {
            result->ndata[row * result->ncol + col] = 0;
            for (size_t k = 0; k < this->ncol; k++) {
                result->ndata[row * result->ncol + col] +=
                        this->ndata[row * this->ncol + k] * mat.ndata[k * mat.ncol + col];
            }
        }
    }
    return *result;
}

Matrix Matrix::operator-(const Matrix &mat) {
    if (this->nrow != mat.nrow || this->ncol != mat.ncol || !this->is_valid || !mat.is_valid) return {};

    auto *result = new Matrix(mat.nrow, mat.ncol);
    if (!result->ndata) return {};

    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            result->data()[row * result->cols() + col] =
                    this->data()[row * this->cols() + col] - mat.data()[row * mat.cols() + col];
        }
    }
    result->is_valid = true;
    return *result;
}

Matrix Matrix::operator+(const Matrix &mat) {
    if (this->nrow != mat.nrow || this->ncol != mat.ncol || !this->is_valid || !mat.is_valid) return {};

    auto *result = new Matrix(mat.nrow, mat.ncol);
    if (!result->ndata) return {};

    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            result->data()[row * result->cols() + col] =
                    this->data()[row * this->cols() + col] + mat.data()[row * mat.cols() + col];
        }
    }
    result->is_valid = true;
    return *result;
}

Matrix Matrix::operator*(double value) {
    if (!this->is_valid) return *this;

    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->data()[row * this->cols() + col] *= value;
        }
    }

    return *this;
}

Matrix Matrix::operator/(double value) {
    if (!this->is_valid) return *this;

    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->coeffRef(row, col) /= value;
        }
    }

    return *this;
}

Matrix &Matrix::operator=(const Matrix &mat) {
    if (this == &mat) return *this;
    delete[] this->ndata;
    this->nrow = mat.nrow;
    this->ncol = mat.ncol;
    this->ndata = new double[this->nrow * this->ncol];

    std::memcpy(this->ndata, mat.ndata, this->nrow * this->ncol);

    this->is_valid = true;
    return *this;
}

Matrix &Matrix::operator*=(const Matrix &mat) {
    *this = *this * mat;
    return *this;
}

Matrix &Matrix::operator+=(const Matrix &mat) {
    *this = *this + mat;
    return *this;
}

Matrix &Matrix::operator-=(const Matrix &mat) {
    *this = *this - mat;
    return *this;
}

Matrix &Matrix::operator*=(double value) {
    *this = *this * value;
    return *this;
}

Matrix &Matrix::operator/=(double value) {
    *this = *this / value;
    return *this;
}

bool Matrix::isValid() const {
    return this->is_valid;
}

void Matrix::resize(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) {
        this->nrow = 0;
        this->ncol = 0;
        this->ndata = nullptr;
    } else {
        delete[] this->ndata;
        this->nrow = rows;
        this->ncol = cols;
        this->ndata = new double[this->nrow * this->ncol];
    }
    this->is_valid = false;
}

const double &Matrix::coeffRef(size_t rowIdx, size_t colIdx) const {
    if (rowIdx > this->nrow) throw std::out_of_range("Row index out of range");
    if (colIdx > this->ncol) throw std::out_of_range("Col index out of range");

    return this->ndata[rowIdx * this->ncol + colIdx];
}

double &Matrix::coeffRef(size_t rowIdx, size_t colIdx) {
    if (rowIdx > this->nrow) throw std::out_of_range("Row index out of range");
    if (colIdx > this->ncol) throw std::out_of_range("Col index out of range");

    return this->ndata[rowIdx * this->ncol + colIdx];
}

const double *Matrix::data() const {
    return this->ndata;
}

double *Matrix::data() {
    return this->ndata;
}

size_t Matrix::rows() const {
    return this->nrow;
}

size_t Matrix::cols() const {
    return this->ncol;
}

Matrix &Matrix::setIdentity() {
    if (this->ndata == nullptr) {
        this->is_valid = false;
        return *this;
    }
    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->ndata[row * this->ncol + col] = (row == col);
        }
    }
    this->is_valid = true;
    return *this;
}

Matrix &Matrix::setZero() {
    if (!this->ndata) {
        this->is_valid = false;
        return *this;
    }
    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->ndata[row * this->ncol + col] = 0;
        }
    }
    this->is_valid = true;
    return *this;
}

Matrix &Matrix::setConstants(double value) {
    if (!this->ndata) {
        this->is_valid = false;
        return *this;
    }
    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->ndata[row * this->ncol + col] = value;
        }
    }
    this->is_valid = true;
    return *this;
}

Matrix &Matrix::setIdentity(size_t rows, size_t cols) {
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

Matrix Matrix::transpose() {
    Matrix temp(this->nrow, this->ncol);

    if (!temp.ndata) {
        this->is_valid = false;
        return *this;
    }

    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            temp.ndata[row * temp.ncol + col] = this->ndata[col * this->ncol + row];
        }
    }
    return *this = temp;
}

Matrix Matrix::inverse() {
    if (!this->is_valid || this->ncol != this->nrow) return {};

    double det = this->det();

    *this = this->transpose() / det;

    return *this;
}

double determinant(double **matrix, size_t size) {
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

double Matrix::det() {
    if (!this->is_valid || this->ncol != this->nrow) return std::nan("NaN");
    size_t size = this->ncol;

    auto **temp = new double *[size];
    for (size_t row = 0; row < size; row++) {
        temp[row] = new double[size];
        for (size_t col = 0; col < size; col++) {
            temp[row][col] = this->ndata[row * this->ncol + col];
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
    return Matrix().setIdentity(rows, cols);
}

Matrix Matrix::zeros(size_t rows, size_t cols) {
    return Matrix().setZero(rows, cols);
}

Matrix Matrix::constants(size_t rows, size_t cols, double value) {
    return Matrix(rows, cols).setConstants(value);
}

Matrix operator*(double value, const Matrix &mat) {
    Matrix temp = mat;
    return temp * value;
}

std::ostream &operator<<(std::ostream &stream, const Matrix &matrix) {
    for (size_t row = 0; row < matrix.nrow; row++) {
        for (size_t col = 0; col < matrix.ncol; col++) {
            stream << "|" << matrix.ndata[row * matrix.ncol + col];
        }
        stream << "|" << std::endl;
    }
    return stream;
}
