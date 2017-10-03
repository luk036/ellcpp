cmake_minimum_required (VERSION 2.8.11)
project (ellcpp)
set ( CMAKE_BUILD_TYPE Debug )
add_definitions ( -Wall -std=c++1z -fcoroutines-ts )
add_executable (Main profit_main.cpp)
target_link_libraries (Main)


