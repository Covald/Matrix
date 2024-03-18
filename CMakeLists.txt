cmake_minimum_required(VERSION 3.22)
project(example_project)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -Wpedantic)

include_directories(src)

add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
#add_subdirectory(lib/googletest)
#add_executable(${CMAKE_PROJECT_NAME}_run main.cpp)
#
## set this flag when running coverage tests in Clion
##set(CMAKE_CXX_FLAGS "--coverage")
#
## generate debug information
#set(CMAKE_CXX_FLAGS "-g")
#
#include_directories(src)
#add_subdirectory(src)
#
#add_subdirectory(tests)
#
#target_link_libraries(${CMAKE_PROJECT_NAME}_run ${CMAKE_PROJECT_NAME}_lib)
