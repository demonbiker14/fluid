cmake_minimum_required(VERSION 3.30)
project(fluid)

set(CMAKE_CXX_STANDARD 20)

add_executable(original fluid.cpp)
add_executable(task2 main.cpp)

