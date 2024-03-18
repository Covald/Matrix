#include "gtest/gtest.h"
#include "../src/ExampleClass.cpp"

TEST(ExampleSuite, ExampleTest) {
    ExampleClass e = ExampleClass();

    EXPECT_EQ(1, e.doSomething());
}