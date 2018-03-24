cmake_minimum_required (VERSION 3.3)

project (ellcpp)

set_property (GLOBAL PROPERTY USE_FOLDERS ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)
add_definitions ( -std=c++17 )
add_subdirectory (third_party EXCLUDE_FROM_ALL)

add_executable (Main profit_main.cpp)
target_link_libraries ( Main -lfmt)

enable_testing ()


