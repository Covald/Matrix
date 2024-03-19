//
// Created by Covald on 19.03.2024.
//
#include <iostream>
#include "Matrix.hpp"

Matrix &set_matrix_by_vector(Matrix &mtx, const double *array, size_t size) {
    if (size != mtx.rows() * mtx.cols()) return mtx;

    for (size_t row = 0; row < mtx.rows(); row++) {
        for (size_t col = 0; col < mtx.cols(); col++) {
            mtx.coeffRef(row, col) = array[row * mtx.cols() + col];
        }
    }

    return mtx;
}

int main(int argc, char **argv) {
    size_t size;
    std::cin >> size;
    char *result = new char[size];
    int a = 'a';
    result[1] = (char) 83;
    result[2] = (char) a;
    std::cout << result[1] << result[2];
    delete[] result;
//    Matrix m(3, 3);
//    double a[] = {10,2,3,4,5,6,7,8,9};
//    m = set_matrix_by_vector(m,a,9);
//    std::cout << m << std::endl;
//    std::cout << m.transpose() << std::endl;
}