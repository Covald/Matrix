#include <cstdlib>
#include <iostream>

class Matrix final {
private:
    size_t nrow, ncol;
    bool is_valid;
    double *ndata;

public:
    Matrix();

    [[maybe_unused]] explicit Matrix(size_t cols);

    Matrix(size_t rows, size_t cols);

    ~Matrix();

    [[maybe_unused]] Matrix(const Matrix &mat);

    Matrix operator*(const Matrix &mat);

    Matrix operator-(const Matrix &mat);

    Matrix operator+(const Matrix &mat);

    Matrix operator*(double value);

    Matrix operator/(double value);

    Matrix &operator=(const Matrix &mat);

    Matrix &operator*=(const Matrix &mat);

    Matrix &operator+=(const Matrix &mat);

    Matrix &operator-=(const Matrix &mat);

    Matrix &operator*=(double value);

    Matrix &operator/=(double value);

    [[nodiscard]] bool isValid() const;

    void resize(size_t rows, size_t cols);

    [[nodiscard]] const double &coeffRef(size_t rowIdx, size_t colIdx) const;

    double &coeffRef(size_t rowIdx, size_t colIdx);

    [[nodiscard]] const double *data() const;

    double *data();

    [[nodiscard]] size_t rows() const;

    [[nodiscard]] size_t cols() const;

    Matrix &setIdentity();

    Matrix &setZero();

    Matrix &setConstants(double value);

    Matrix &setIdentity(size_t rows, size_t cols);

    Matrix &setZero(size_t rows, size_t cols);

    Matrix &setConstants(size_t rows, size_t cols, double value);

    Matrix transpose();

    Matrix inverse();

    double det();

    static Matrix identity(size_t rows, size_t cols);

    static Matrix zeros(size_t rows, size_t cols);

    static Matrix constants(size_t rows, size_t cols, double value);

    friend Matrix operator*(double value, const Matrix &mat);

    friend std::ostream &operator<<(std::ostream &stream, const Matrix &matrix);
};
