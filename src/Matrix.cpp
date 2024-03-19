//
// Created by Covald on 19.03.2024.
//

#include "Matrix.hpp"
#include <stdexcept>


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
    this->nrow = rows;
    this->ncol = cols;
    this->ndata = new double[this->nrow * this->ncol];
    this->is_valid = false;
}

Matrix::~Matrix() {
    delete[] this->ndata;
}

[[maybe_unused]] Matrix::Matrix(const Matrix &mat) {
    this->nrow = mat.nrow;
    this->ncol = mat.ncol;
    this->ndata = new double[this->nrow * this->ncol];

    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->coeffRef(row, col) = mat.coeffRef(row, col);
        }
    }

    this->is_valid = true;
}

Matrix Matrix::operator*(const Matrix &mat) {
    if (this->nrow != mat.ncol || !this->is_valid || !mat.is_valid) return {};

    auto *result = new Matrix(this->nrow, mat.ncol);
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

    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->coeffRef(row, col) = mat.coeffRef(row, col);
        }
    }

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
    if (cols < 1 || rows < 1) throw std::invalid_argument("Arguments must be greater than 0.");

    delete[] this->ndata;
    this->nrow = rows;
    this->ncol = cols;
    this->ndata = new double[this->nrow * this->ncol];
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
    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->ndata[row * this->ncol + col] = (row == col);
        }
    }
    this->is_valid = true;
    return *this;
}

Matrix &Matrix::setZero() {
    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            this->ndata[row * this->ncol + col] = 0;
        }
    }
    this->is_valid = true;
    return *this;
}

Matrix &Matrix::setConstants(double value) {
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
    for (size_t row = 0; row < this->nrow; row++) {
        for (size_t col = 0; col < this->ncol; col++) {
            temp.ndata[row * temp.ncol + col] = this->ndata[col * this->ncol + row];
        }
    }
    return *this = temp;
}

Matrix Matrix::inverse() {
    return Matrix();
}

double Matrix::det() {
    return 0;
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
    return Matrix();
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
