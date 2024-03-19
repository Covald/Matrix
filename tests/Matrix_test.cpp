#include "gtest/gtest.h"
#include "../src/Matrix.cpp"
#include "stdexcept"
#include <iostream>

TEST(Contructors, base_contructor) {
    Matrix m;
    ASSERT_EQ(m.isValid(), false);
}