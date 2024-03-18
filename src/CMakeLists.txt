set(BINARY ${CMAKE_PROJECT_NAME})

#add_executable(ExampleProject main.cpp example.h example.cpp)
file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.hpp *.cpp)
set(SOURCES ${SOURCES})

add_executable(${BINARY}_run ${SOURCES})

add_library(${BINARY}_lib STATIC ${SOURCES})