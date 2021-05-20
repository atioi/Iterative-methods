cmake_minimum_required(VERSION 3.19)
project(laboratorium7)

set(CMAKE_CXX_STANDARD 14)

add_executable(laboratorium7 main.cpp estimator.cpp methods.h other.h SOR.cpp other.cpp estimator.h residuum.cpp residuum.h Jacobi.cpp Gauss_Seidel.cpp )