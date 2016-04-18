cmake_minimum_required(VERSION 3.0)
find_package(Eigen3 REQUIRED)

set (CMAKE_CXX_STANDARD 11)

include_directories(${EIGEN3_INCLUDE_DIR})


add_library(TimexDSpline TimexDSpline.cpp)
add_executable(spline_test spline_test.cpp) 
target_link_libraries(spline_test TimexDSpline)
