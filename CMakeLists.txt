# Distributed under the MIT License (See accompanying file /LICENSE )

# CMake build : global project

cmake_minimum_required (VERSION 3.3)

project (Ell)

set_property (GLOBAL PROPERTY USE_FOLDERS ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set(LIBS ${LIBS} "-L$ENV{CONDA_PREFIX}/lib")

set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)
find_package (xtensor REQUIRED)
message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}")
# set(LIBS ${LIBS} "-L${xtensor_INCLUDE_DIRS}/../lib")
# message(STATUS "Test: ${LIBS}")

find_package (xtensor-blas REQUIRED)
message(STATUS "Found xtensor-blas: ${xtensor_blas_INCLUDE_DIRS}")
# set(LIBS ${LIBS} "-L${xtensor_INCLUDE_DIRS}/../lib")
# message(STATUS "Test: ${LIBS}")

# set (OpenBLAS_DIR $ENV{CONDA}/lib/cmake/openblas)
find_package (OpenBLAS REQUIRED)
message(STATUS "Found openblas: ${OpenBLAS_LIBRARIES}")
set(LIBS ${LIBS} ${OpenBLAS_LIBRARIES})
message(STATUS "LIB: ${LIBS}")

# find_package (lapack REQUIRED)
# message(STATUS "Found lapack: ${lapack_INCLUDE_DIRS}")
# set(LIBS ${LIBS} ${lapack_LIBRARIES})
# message(STATUS "LIB: ${LIBS}")

add_definitions ( -std=c++1z -g)
add_subdirectory (third_party EXCLUDE_FROM_ALL)

# add_executable (Main profit_main.cpp)
# add_executable (Test exp1d.cpp)
# target_link_libraries ( Main Threads::Threads ${LIBS} -llapack -lblas -lfmt)
# target_link_libraries ( Test Threads::Threads ${LIBS} -llapack -lblas -lfmt)
set (LIBRARY_INCLUDE_PATH  ${LIBRARY_INCLUDE_PATH} ${xtensor_INCLUDE_DIRS})

add_subdirectory (lib)
add_subdirectory (app)

enable_testing ()
