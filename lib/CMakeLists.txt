# Distributed under the MIT License (See accompanying file /LICENSE )

# Adds convenience methods, see cmake/cleanCppExtensions.cmake
include(cleanCppExtensions)

# Creates options to turn on sanitizers, see cmake/sanitizers.cmake
include(sanitizers)

# CMake build : library

#configure variables
set (LIB_NAME "${PROJECT_NAME}Lib")

#configure directories
set (LIBRARY_MODULE_PATH "${PROJECT_SOURCE_DIR}/lib")
set (LIBRARY_SRC_PATH  "${LIBRARY_MODULE_PATH}/src" )
set (LIBRARY_INCLUDE_PATH  "${LIBRARY_MODULE_PATH}/include")

#set includes
# include_directories (${LIBRARY_INCLUDE_PATH} ${THIRD_PARTY_INCLUDE_PATH})
include_directories (${LIBRARY_INCLUDE_PATH} ${THIRD_PARTY_INCLUDE_PATH} ${xtensor_INCLUDE_DIRS})

#set sources
file (GLOB LIB_HEADER_FILES "${LIBRARY_INCLUDE_PATH}/*.hpp ${LIBRARY_INCLUDE_PATH}/*.h")
file (GLOB LIB_SOURCE_FILES "${LIBRARY_SRC_PATH}/*.cpp")

#set library
add_library (${LIB_NAME} STATIC ${LIB_SOURCE_FILES} ${LIB_HEADER_FILES})

#export vars
set (LIBRARY_INCLUDE_PATH  ${LIBRARY_INCLUDE_PATH} PARENT_SCOPE)
set (LIB_NAME ${LIB_NAME} PARENT_SCOPE)

#test
enable_testing ()
add_subdirectory (test)

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if( NOT cmake_build_type_tolower STREQUAL "coverage"
    AND cmake_build_type_tolower STREQUAL "release")
    # If you want to self-test benchmark lib too, turn me ON
    #
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable benchmark self-testing?")
    
    # Fix dependency download
    # https://github.com/google/benchmark/issues/564
    #
    set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON CACHE BOOL "Download dependencies?")
    include(ConfigGBench)

    add_subdirectory (bench)
endif()


#------------------------------------------------------------------------------
# Other MISC targets

# Adds misc targets: format, cppcheck, tidy, see cmake/cleanCppExtensions.cmake
addMiscTargets()
