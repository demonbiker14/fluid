cmake_minimum_required(VERSION 3.30)
project(fluid)

set(CMAKE_CXX_STANDARD 20)

add_executable(
        original fluid.cpp
)

set(DEFAULT_TYPES FLOAT,DOUBLE,FIXED(32,16),FAST_FIXED(16,8))
set(TYPES ${DEFAULT_TYPES} CACHE STRING "Specify the TYPES for the simulation")

add_executable(task2 main.cpp)
target_compile_definitions(task2 PRIVATE TYPES=${TYPES})
add_executable(task3 main.cpp)

