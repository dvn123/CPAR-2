cmake_minimum_required(VERSION 3.1)
project(CPAR_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(SOURCE_FILES blocksOpenMP.cpp)
add_executable(CPAR_2 ${SOURCE_FILES})