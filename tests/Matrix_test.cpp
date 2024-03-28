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

TEST(Constructors, default_constructor) {
    Matrix m;
    ASSERT_EQ(m.rows(), 0);
    ASSERT_EQ(m.cols(), 0);
    ASSERT_FALSE(m.data());
    ASSERT_FALSE(m.isValid());
}

TEST(Constructors, two_arguments) {
    Matrix m(3, 3);
    std::cout << m << std::endl;
    ASSERT_EQ(m.rows(), 3);
    ASSERT_EQ(m.cols(), 3);
    ASSERT_TRUE(m.data());
    ASSERT_TRUE(m.isValid());
    for (size_t row = 0; row < 3; row++)
        for (size_t col = 0; col < 3; col++)
            ASSERT_EQ(m.coeffRef(row, col), 0);
}

TEST(Constructors, two_bad_arguments) {
    Matrix m(0, 0);
    std::cout << m << std::endl;
    ASSERT_EQ(m.rows(), 0);
    ASSERT_EQ(m.cols(), 0);
    ASSERT_FALSE(m.data());
    ASSERT_FALSE(m.isValid());
}

TEST(Constructors, one_argument) {
    Matrix m(3);
    std::cout << m << std::endl;
    ASSERT_EQ(m.rows(), 1);
    ASSERT_EQ(m.cols(), 3);
    ASSERT_TRUE(m.data());
    ASSERT_TRUE(m.isValid());
    for (size_t col = 0; col < 3; col++) {
        ASSERT_EQ(m.coeffRef(0, col), 0);
    }
}

TEST(Constructors, one_bad_argument) {
    Matrix m(0);
    std::cout << m << std::endl;
    ASSERT_EQ(m.rows(), 0);
    ASSERT_EQ(m.cols(), 0);
    ASSERT_FALSE(m.data());
    ASSERT_FALSE(m.isValid());
}

TEST(Constructors, copy) {
    Matrix m1(3, 3);
    double array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    set_matrix_by_vector(m1, array, 9);
    std::cout << m1 << std::endl;
    for (size_t row = 0; row < 3; row++)
        for (size_t col = 0; col < 3; col++)
            ASSERT_EQ(m1.coeffRef(row, col), array[row * 3 + col]);

    Matrix m2 = m1;
    std::cout << m2 << std::endl;
    ASSERT_EQ(m2.rows(), m1.rows());
    ASSERT_EQ(m2.cols(), m1.cols());
    ASSERT_TRUE(m2.data());
    ASSERT_TRUE(m2.isValid());
    for (size_t row = 0; row < 3; row++)
        for (size_t col = 0; col < 3; col++)
            ASSERT_EQ(m2.coeffRef(row, col), m1.coeffRef(row, col));
}

TEST(BaseFunctions, resize) {
    Matrix m(3, 3);
    m.resize(1, 1);
    ASSERT_EQ(m.rows(), 1);
    ASSERT_EQ(m.cols(), 1);
    ASSERT_TRUE(m.isValid());
    ASSERT_TRUE(m.data());
    ASSERT_EQ(m.coeffRef(0, 0), 0);
}

TEST(BaseFunctions, coeff_ref) {
    Matrix m(3, 3);
    m.coeffRef(0, 0) = 1;
    ASSERT_EQ(m.coeffRef(0, 0), 1);
}

TEST(BaseFunctions, bad_coeff_ref) {
    Matrix m(3, 3);
    ASSERT_THROW(m.coeffRef(3, 3), std::out_of_range);
}

TEST(BaseFunctions, set_identity) {
    Matrix m = Matrix::identity(3, 3);
    std::cout << m << std::endl;
    for (size_t row = 0; row < 3; row++)
        for (size_t col = 0; col < 3; col++)
            if (row == col) ASSERT_EQ(m.coeffRef(row, col), 1);
            else
                ASSERT_EQ(m.coeffRef(row, col), 0);
}

TEST(BaseFunctions, set_constants) {
    Matrix m = Matrix::constants(3, 3, 3);
    std::cout << m << std::endl;
    for (size_t row = 0; row < 3; row++)
        for (size_t col = 0; col < 3; col++)
            ASSERT_EQ(m.coeffRef(row, col), 3);
}

TEST(Operators, assigment) {
    Matrix m1(3, 3);
    double array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    set_matrix_by_vector(m1, array, 9);
    std::cout << m1 << std::endl;
    for (size_t row = 0; row < 3; row++)
        for (size_t col = 0; col < 3; col++)
            ASSERT_EQ(m1.coeffRef(row, col), array[row * 3 + col]);

    Matrix m2(3, 3);
    m2 = m1;

    std::cout << m2 << std::endl;
    ASSERT_EQ(m2.rows(), m1.rows());
    ASSERT_EQ(m2.cols(), m1.cols());
    ASSERT_TRUE(m2.data());
    ASSERT_TRUE(m2.isValid());
    for (size_t row = 0; row < 3; row++)
        for (size_t col = 0; col < 3; col++)
            ASSERT_EQ(m2.coeffRef(row, col), m1.coeffRef(row, col));
}

TEST(Operators, multiply_mat) {
    Matrix m1(3, 2);
    Matrix m2(2, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6};
    double m2_array[] = {1, 2, 3, 4, 5, 6};
    set_matrix_by_vector(m1, m1_array, 6);
    set_matrix_by_vector(m2, m2_array, 6);
    Matrix m3 = m1 * m2;
    std::cout << m1 << std::endl
              << m2 << std::endl
              << m3 << std::endl;
    ASSERT_EQ(m3.rows(), m1.rows());
    ASSERT_EQ(m3.cols(), m2.cols());
    ASSERT_TRUE(m3.isValid());
    double m3_array[] = {9, 12, 15, 19, 26, 33, 29, 40, 51};
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m2.cols(); col++)
            ASSERT_EQ(m3.coeffRef(row, col), m3_array[row * m2.cols() + col]);
}

TEST(Operators, multiply_mat_with_assigment) {
    Matrix m1(3, 2);
    Matrix m2(2, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6};
    double m2_array[] = {1, 2, 3, 4, 5, 6};
    set_matrix_by_vector(m1, m1_array, 6);
    set_matrix_by_vector(m2, m2_array, 6);
    std::cout << m1 << std::endl
              << m2 << std::endl;
    m1 *= m2;
    std::cout << m1 << std::endl;
    ASSERT_EQ(m1.rows(), m1.rows());
    ASSERT_EQ(m1.cols(), m2.cols());
    ASSERT_TRUE(m1.isValid());
    double m3_array[] = {9, 12, 15, 19, 26, 33, 29, 40, 51};
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m2.cols(); col++)
            ASSERT_EQ(m1.coeffRef(row, col), m3_array[row * m2.cols() + col]);
}

TEST(Operators, addition_mat) {
    Matrix m1(3, 3);
    Matrix m2(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m2_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m3_array[] = {2, 4, 6, 8, 10, 12, 14, 16, 18};
    set_matrix_by_vector(m1, m1_array, 9);
    set_matrix_by_vector(m2, m2_array, 9);
    Matrix m3 = m1 + m2;
    std::cout << m1 << m2 << m3 << std::endl;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m2.cols(); col++)
            ASSERT_EQ(m3.coeffRef(row, col), m3_array[row * m2.cols() + col]);
}

TEST(Operators, bad_addition_mat) {
    Matrix m1(3, 2);
    Matrix m2(2, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6};
    double m2_array[] = {1, 2, 3, 4, 5, 6};
    set_matrix_by_vector(m1, m1_array, 6);
    set_matrix_by_vector(m2, m2_array, 6);
    Matrix m3 = m1 + m2;
    std::cout << m1 << m2 << m3 << std::endl;
    ASSERT_FALSE(m3.isValid());
}

TEST(Operators, addition_mat_with_assigment) {
    Matrix m1(3, 3);
    Matrix m2(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m2_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m3_array[] = {2, 4, 6, 8, 10, 12, 14, 16, 18};
    set_matrix_by_vector(m1, m1_array, 9);
    set_matrix_by_vector(m2, m2_array, 9);
    std::cout << m1 << m2 << std::endl;
    m1 += m2;
    std::cout << m1 << std::endl;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m2.cols(); col++)
            ASSERT_EQ(m1.coeffRef(row, col), m3_array[row * m2.cols() + col]);
}

TEST(Operators, bad_addition_mat_with_assigment) {
    Matrix m1(3, 2);
    Matrix m2(2, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6};
    double m2_array[] = {1, 2, 3, 4, 5, 6};
    set_matrix_by_vector(m1, m1_array, 6);
    set_matrix_by_vector(m2, m2_array, 6);
    std::cout << m1 << m2 << std::endl;
    m1 += m2;
    std::cout << m1 << std::endl;
    ASSERT_FALSE(m1.isValid());
}

TEST(Operators, subtraction_mat) {
    Matrix m1(3, 3);
    Matrix m2(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m2_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m3_array[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    set_matrix_by_vector(m1, m1_array, 9);
    set_matrix_by_vector(m2, m2_array, 9);
    Matrix m3 = m1 - m2;
    std::cout << m1 << m2 << m3 << std::endl;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m2.cols(); col++)
            ASSERT_EQ(m3.coeffRef(row, col), m3_array[row * m2.cols() + col]);
}

TEST(Operators, bad_subtraction_mat) {
    Matrix m1(3, 2);
    Matrix m2(2, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6};
    double m2_array[] = {1, 2, 3, 4, 5, 6};
    set_matrix_by_vector(m1, m1_array, 6);
    set_matrix_by_vector(m2, m2_array, 6);
    Matrix m3 = m1 - m2;
    std::cout << m1 << m2 << m3 << std::endl;
    ASSERT_FALSE(m3.isValid());
}

TEST(Operators, subtraction_mat_with_assigment) {
    Matrix m1(3, 3);
    Matrix m2(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m2_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double m3_array[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    set_matrix_by_vector(m1, m1_array, 9);
    set_matrix_by_vector(m2, m2_array, 9);
    std::cout << m1 << m2 << std::endl;
    m1 -= m2;
    std::cout << m1 << std::endl;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m2.cols(); col++)
            ASSERT_EQ(m1.coeffRef(row, col), m3_array[row * m2.cols() + col]);
}

TEST(Operators, bad_subtraction_mat_with_assigment) {
    Matrix m1(3, 2);
    Matrix m2(2, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6};
    double m2_array[] = {1, 2, 3, 4, 5, 6};
    set_matrix_by_vector(m1, m1_array, 6);
    set_matrix_by_vector(m2, m2_array, 6);
    std::cout << m1 << m2 << std::endl;
    m1 -= m2;
    std::cout << m1 << std::endl;
    ASSERT_FALSE(m1.isValid());
}

TEST(Operators_numeric, multiply) {
    Matrix m1(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    set_matrix_by_vector(m1, m1_array, 9);
    double k = 3;
    double result[] = {3, 6, 9, 12, 15, 18, 21, 24, 27};
    Matrix m2 = m1 * k;
    std::cout << m1 << k << std::endl << m2;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m1.cols(); col++)
            ASSERT_EQ(m2.coeffRef(row, col), result[row * m1.cols() + col]);
}

TEST(Operators_numeric, divide) {
    Matrix m1(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double k = 3;
    double result[] = {3, 6, 9, 12, 15, 18, 21, 24, 27};
    set_matrix_by_vector(m1, result, 9);

    Matrix m2 = m1 / k;
    std::cout << m1 << k << std::endl << m2;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m1.cols(); col++)
            ASSERT_EQ(m2.coeffRef(row, col), m1_array[row * m1.cols() + col]);
}

TEST(Operators_numeric, multiply_with_assigment) {
    Matrix m1(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    set_matrix_by_vector(m1, m1_array, 9);
    double k = 3;
    double result[] = {3, 6, 9, 12, 15, 18, 21, 24, 27};

    std::cout << m1 << k << std::endl;
    m1 *= k;
    std::cout << m1;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m1.cols(); col++)
            ASSERT_EQ(m1.coeffRef(row, col), result[row * m1.cols() + col]);
}

TEST(Operators_numeric, divide_with_assigment) {
    Matrix m1(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double k = 3;
    double result[] = {3, 6, 9, 12, 15, 18, 21, 24, 27};
    set_matrix_by_vector(m1, result, 9);
    std::cout << m1 << k << std::endl;
    m1 /= k;
    std::cout << m1;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m1.cols(); col++)
            ASSERT_EQ(m1.coeffRef(row, col), m1_array[row * m1.cols() + col]);
}

TEST(Complex_functions, transpose) {
    Matrix m1(3, 3);
    double m1_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    set_matrix_by_vector(m1, m1_array, 9);
    Matrix m2 = m1.transpose();
    std::cout << m1 << m2 << std::endl;
    for (size_t row = 0; row < m1.rows(); row++)
        for (size_t col = 0; col < m1.cols(); col++)
            ASSERT_EQ(m2.coeffRef(col, row), m1_array[row * m1.cols() + col]);
}

TEST(Complex_functions, inverse) {
    Matrix m1(3, 3);
    double m1_array[] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
    double m2_array[] = {1, -1, 1, -38, 41, -34, 27, -29, 24};
    set_matrix_by_vector(m1, m1_array, 9);
    Matrix m2 = m1.inverse();
    std::cout << m1 << std::endl << m2 << std::endl;
    for (size_t row = 0; row < m2.rows(); row++)
        for (size_t col = 0; col < m2.cols(); col++)
            EXPECT_LE(
                    std::abs(m2_array[row * m2.cols() + col] - m2.coeffRef(row, col)),
                    1E-06);
}

TEST(Complex_functions, determinant) {
    Matrix m1(3, 3);
    double m1_array[] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
    set_matrix_by_vector(m1, m1_array, 9);
    std::cout << m1 << std::endl;
    double det = m1.det();
    std::cout << "determinant=" << det << std::endl;
    ASSERT_EQ(det, -1);
}