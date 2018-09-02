cmake_minimum_required(VERSION 3.12)
project(Mandlebrot_MPI)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_CXX_STANDARD 14)

add_executable(Mandlebrot_MPI main.cpp)