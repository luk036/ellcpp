cmake_minimum_required (VERSION 3.3)

project (ellcpp)

set_property (GLOBAL PROPERTY USE_FOLDERS ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)
find_package (xtensor REQUIRED)
message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}")
set(LIBS ${LIBS} "-L${xtensor_INCLUDE_DIRS}/../lib")
message(STATUS "Test: ${LIBS}")

set (LIBRARY_INCLUDE_PATH  ${LIBRARY_INCLUDE_PATH} ${xtensor_INCLUDE_DIRS})

add_definitions ( -std=c++17 )
add_subdirectory (third_party EXCLUDE_FROM_ALL)

include_directories (${LIBRARY_INCLUDE_PATH} ${TEST_THIRD_PARTY_INCLUDE_PATH})
add_executable (Main profit_main.cpp)
add_executable (Test exp1d.cpp)
target_link_libraries ( Main Threads::Threads ${LIBS} -llapack -lblas -lfmt)
target_link_libraries ( Test Threads::Threads ${LIBS} -llapack -lblas -lfmt)

enable_testing ()


