# Distributed under the MIT License (See accompanying file /LICENSE )

# CMake build : global project

cmake_minimum_required (VERSION 3.3)

project (Ell)

set_property (GLOBAL PROPERTY USE_FOLDERS ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# set(LIBS ${LIBS} "-L$ENV{CONDA_PREFIX}/lib")

set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)
find_package (xtensor-blas REQUIRED)
if (xtensor-blas_FOUND)
    message(STATUS "Found xtensor-blas: ${xtensor-blas_INCLUDE_DIRS}")
    # set(LIBS ${LIBS} ${xtensor-blas_LIBRARIES})
endif (xtensor-blas_FOUND)
find_package (xtensor REQUIRED)
if (xtensor_FOUND)
    message(STATUS "Found xtensor: ${xtensor_INCLUDE_DIRS}")
    set(LIBS ${LIBS} ${xtensor_LIBRARIES})
endif (xtensor_FOUND)
find_package (OpenBLAS REQUIRED)
if (OpenBLAS_FOUND)
    message(STATUS "Found OpenBLAS: ${OpenBLAS_INCLUDE_DIRS}")
    set(LIBS ${LIBS} ${OpenBLAS_LIBRARIES})
endif (OpenBLAS_FOUND)
find_package (Boost REQUIRED)
if (Boost_FOUND)
    message(STATUS "Found boost: ${Boost_INCLUDE_DIRS}")
    # set(LIBS ${LIBS} ${Boost_LIBRARIES})
endif (Boost_FOUND)

# add_definitions ( -std=c++1z -g)
add_subdirectory (third_party EXCLUDE_FROM_ALL)

# add_executable (Main profit_main.cpp)
# add_executable (Test exp1d.cpp)
# target_link_libraries ( Main Threads::Threads ${LIBS} -llapack -lblas -lfmt)
# target_link_libraries ( Test Threads::Threads ${LIBS} -llapack -lblas -lfmt)
set (LIBRARY_INCLUDE_PATH  ${LIBRARY_INCLUDE_PATH} ${xtensor_INCLUDE_DIRS})

add_subdirectory (lib)
add_subdirectory (app)

enable_testing ()
