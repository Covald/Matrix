#include "gtest/gtest.h"
#include "../src/Matrix.cpp"
#include "stdexcept"
#include <iostream>

Matrix &set_matrix_by_vector(Matrix &mtx, const double *array, size_t size) {
    if (size != mtx.rows() * mtx.cols()) return mtx;

    for (size_t row = 0; row < mtx.rows(); row++) {
        for (size_t col = 0; col < mtx.cols(); col++) {
            mtx.coeffRef(row, col) = array[row * mtx.cols() + col];
        }
    }

    return mtx;
}

TEST(Contructors, base_contructor) {
    Matrix m;
    ASSERT_EQ(m.isValid(), false);
}