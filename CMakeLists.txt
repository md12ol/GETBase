cmake_minimum_required(VERSION 3.16)
project(THADSNBit)

set(CMAKE_CXX_STANDARD 20)

add_executable(THADSNBit
        main.cpp
        setu.h
        setu.cpp
        stat.h
        stat.cpp
        bitspray.cpp
        bitspray.h)

set_property(TARGET THADSNBit PROPERTY WORKING_DIRECTORY ./)