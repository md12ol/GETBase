cmake_minimum_required(VERSION 3.16)
project(GETBase)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_BUILD_TYPE "Release")

#add_executable(GETBase
#        main.cpp
#        main.h
#        SDA.cpp
#        SDA.h
#        Graph.h
#        Graph.cpp)

#add_executable(GETBase
#       GraphTest.cpp
#        Graph.cpp
#        Graph.h)

add_executable(GETBase
SDATest.cpp
SDA.cpp
SDA.h)

set_property(TARGET GETBase PROPERTY WORKING_DIRECTORY ./)