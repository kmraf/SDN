cmake_minimum_required(VERSION 3.21)
project(graph)

set(CMAKE_CXX_STANDARD 14)

add_executable(graph
        dist_func.cpp
        main.cpp
        other_algs.cpp
        parse_and_wright.cpp)
