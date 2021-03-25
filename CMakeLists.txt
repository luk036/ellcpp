# Distributed under the MIT License (See accompanying file /LICENSE )

# CMake build : global project

# Project specific options :
#   - BP_USE_DOXYGEN
#   - BP_BUILD_TESTS (requires BUILD_TESTING set to ON)
# Other options might be available through the cmake scripts including (not exhaustive):
#   - ENABLE_WARNINGS_SETTINGS
#   - ENABLE_LTO
#
cmake_minimum_required(VERSION 3.8.2)

set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX})

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Do not build in-source. Please remove CMakeCache.txt and the CMakeFiles/ directory. Then build out-of-source.")
endif()

project (EllCpp LANGUAGES CXX)

option(ENABLE_COVERAGE "Generates the coverage build" OFF)
option(ENABLE_DOCTESTS "Include tests in the library. Setting this to OFF will remove all doctest related code.
                        Tests in tests/*.cpp will still be enabled." OFF)

include(CTest) # Must be called before adding tests but after calling project(). This automatically calls enable_testing() and configures ctest targets when using Make/Ninja
include(CMakeDependentOption)# This is a really useful scripts that creates options that depends on other options. It can even be used with generator expressions !

# Custom modules and scripts

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(Colors)
include(Doctest)
include(LTO)
include(Warnings)
include(CopyDllsForDebug)

# Guard against in-source builds and bad build-type strings
#
include(ConfigSafeGuards)
include(CMakeFindDependencyMacro) # seems for xtensor-blas only

###############
##  OPTIONS  ##
###############

# You should try to give as much control over the project setup to the user.
# When modifying compile flags for example, if they are not mandatory, provide an option.

# Use your own option for tests, in case people use your library through add_subdirectory
cmake_dependent_option(BP_BUILD_TESTS
    "Enable Boilerplate project tests targets" ON # By default we want tests if CTest is enabled
    "BUILD_TESTING" OFF # Stay coherent with CTest variables
)

# It is always easier to navigate in an IDE when projects are organized in folders.
set_property (GLOBAL PROPERTY USE_FOLDERS ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)
if(MSVC)
	add_compile_options(/std:c++latest /Wv:18 /W3 /WX /wd4819 /wd4146 /D_SILENCE_CXX17_OLD_ALLOCATOR_MEMBERS_DEPRECATION_WARNING)
else()
    # add_compile_options(-fno-exceptions -DXTENSOR_DISABLE_EXCEPTIONS -DBOOST_NO_EXCEPTION -Wall -Wnon-virtual-dtor -Wunused -pedantic)
    add_compile_options(-Wall -Wnon-virtual-dtor -Wunused -pedantic)
endif()

# set(LIBS ${LIBS} "-L$ENV{CONDA_PREFIX}/lib")

set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)
find_package (fmt REQUIRED)
# find_package (spdlog REQUIRED)

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

find_package (xtensor-fftw REQUIRED)
if (xtensor-fftw_FOUND)
    message(STATUS "Found xtensor-fftw: ${xtensor-fftw_INCLUDE_DIRS}")
    # set(LIBS ${LIBS} ${xtensor-blas_LIBRARIES})
endif (xtensor-fftw_FOUND)

# .. fftw
if(MSVC)
    # no long double component, since in the Windows conda-forge build it is not available
    # and the "official" prebuilt long double library can only be used from MinGW
    find_package(FFTW REQUIRED
            COMPONENTS FLOAT_LIB DOUBLE_LIB)
    add_definitions(-DFFTW_NO_LONGDOUBLE)
else()
    find_package(FFTW REQUIRED
            COMPONENTS FLOAT_LIB DOUBLE_LIB LONGDOUBLE_LIB)
endif(MSVC)
if (FFTW_FOUND)
    set(LIBS ${LIBS} ${FFTW_LIBRARIES})
    message(STATUS "Found FFTW: ${FFTW_LIBRARIES}")
    include_directories(${FFTW_INCLUDE_DIRS})
endif (FFTW_FOUND)

# set(OpenBLAS_USE_STATIC TRUE)
# find_package (OpenBLAS REQUIRED)
# if (OpenBLAS_FOUND)
#     if (NOT MSVC)
#         if (EXISTS ${OpenBLAS_LIBRARIES})
#             message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")
#         else() # workaround the bug in miniconda3
#             set(OpenBLAS_LIBRARIES ${CMAKE_PREFIX_PATH}/lib/libopenblas.so.0)
#             if (EXISTS ${OpenBLAS_LIBRARIES})
#                 message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")
#             endif()
#         endif()
#         set(LIBS ${LIBS} ${OpenBLAS_LIBRARIES})
#     endif()
# endif (OpenBLAS_FOUND)

find_package(OpenBLAS REQUIRED)
if (OpenBLAS_FOUND)
    message(STATUS "Found OpenBLAS: ${OpenBLAS_LIBRARIES}")
    set(LIBS ${LIBS} ${OpenBLAS_LIBRARIES})
endif (OpenBLAS_FOUND)

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
    message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
    set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)

find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
    set(LIBS ${LIBS} ${BLAS_LIBRARIES})
endif (BLAS_FOUND)

find_package (Boost REQUIRED COMPONENTS context coroutine)
if (Boost_FOUND)
    message(STATUS "Found boost: ${Boost_INCLUDE_DIRS}")
    set(LIBS ${LIBS} ${Boost_LIBRARIES})
endif (Boost_FOUND)

# add_definitions ( -std=c++1z -g)
add_subdirectory (external EXCLUDE_FROM_ALL)

# set (LIBRARY_INCLUDE_PATH ${LIBRARY_INCLUDE_PATH} ${xtensor_INCLUDE_DIRS})

# Enable code coverage with -DENABLE_COVERAGE=1
if(ENABLE_COVERAGE)
    set(CMAKE_BUILD_TYPE "Coverage")
    include(CodeCoverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_NAME}_test coverage)
endif()

# add_subdirectory (doc)
add_subdirectory (lib)
add_subdirectory (app)

enable_testing ()
