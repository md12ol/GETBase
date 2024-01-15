cmake_minimum_required(VERSION 3.16)
project(GETBase)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_BUILD_TYPE "Release")

add_executable(GETBase
        main.cpp
        main.h
        SDA.cpp
        SDA.h
        Graph.h
        Graph.cpp)

set_property(TARGET GETBase PROPERTY WORKING_DIRECTORY ./)