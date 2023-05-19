cmake_minimum_required(VERSION 3.16)
project(GETBase)

set(CMAKE_CXX_STANDARD 20)

add_executable(GETBase
        main.cpp
        setu.h
        setu.cpp
        stat.h
        stat.cpp
        SDA.cpp
        SDA.h)

set_property(TARGET GETBase PROPERTY WORKING_DIRECTORY ./)