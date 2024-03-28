#include <cstdlib>
#include <iostream>

class Matrix final {
private:
    size_t n_row{}, n_col{};
    bool is_valid{false};
    double *n_data{};

    double &at(size_t row, size_t col);

    double &at(size_t row, size_t col) const;

    bool init(size_t rows, size_t cols);

public:
    Matrix() = default;

    explicit Matrix(size_t cols);

    Matrix(size_t rows, size_t cols);

    ~Matrix();

    Matrix(const Matrix &mat);

    Matrix operator*(const Matrix &mat) const;

    Matrix operator-(const Matrix &mat) const;

    Matrix operator+(const Matrix &mat) const;

    Matrix operator*(double value) const;

    Matrix operator/(double value) const;

    Matrix &operator=(const Matrix &mat);

    Matrix &operator*=(const Matrix &mat);

    Matrix &operator+=(const Matrix &mat);

    Matrix &operator-=(const Matrix &mat);

    Matrix &operator*=(double value);

    Matrix &operator/=(double value);

    bool isValid() const;

    void resize(size_t rows, size_t cols);

    const double &coeffRef(size_t rowIdx, size_t colIdx) const;

    double &coeffRef(size_t rowIdx, size_t colIdx);

    const double *data() const;

    double *data();

    size_t rows() const;

    size_t cols() const;

    Matrix &setIdentity();

    Matrix &setZero();

    Matrix &setConstants(double value);

    Matrix &setIdentity(size_t rows, size_t cols);

    Matrix &setZero(size_t rows, size_t cols);

    Matrix &setConstants(size_t rows, size_t cols, double value);

    Matrix transpose() const;

    Matrix inverse() const;

    double det() const;

    static Matrix identity(size_t rows, size_t cols);

    static Matrix zeros(size_t rows, size_t cols);

    static Matrix constants(size_t rows, size_t cols, double value);

    friend Matrix operator*(double value, const Matrix &mat);

    Matrix submatrix(uint i, uint j) const;

    double minordet(uint i, uint j) const;

    double cofactor(uint i, uint j) const;

    Matrix minor_matrix() const;

    Matrix cofactor_matrix() const;

    Matrix adj() const;
};

std::ostream &operator<<(std::ostream &stream, const Matrix &matrix);