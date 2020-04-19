# Distributed under the MIT License (See accompanying file /LICENSE )

# CMake build : third party

#configure directories
set (THIRD_PARTY_MODULE_PATH "${PROJECT_SOURCE_DIR}/external")

# -------------------------------------------------------------------------
# fmt

#configure directories
set (GSL_MODULE_PATH "${THIRD_PARTY_MODULE_PATH}/GSL")
set (GSL_INCLUDE_PATH "${GSL_MODULE_PATH}/include")

#include custom cmake function
## include ( "${CATCH_MODULE_PATH}/contrib/ParseAndAddCatchTests.cmake")

# -------------------------------------------------------------------------
# benchmark

#configure directories
# set (BENCHMARK_MODULE_PATH "${THIRD_PARTY_MODULE_PATH}/benchmark")
# set (BENCHMARK_INCLUDE_PATH "${BENCHMARK_MODULE_PATH}/include")
# set (LIBS ${LIBS} ${BENCHMARK_MODULE_PATH}/build/src/libbenchmark.a)

# -------------------------------------------------------------------------

#set variables
set (THIRD_PARTY_INCLUDE_PATH  ${GSL_INCLUDE_PATH})

#set variables for tests
set (TEST_THIRD_PARTY_INCLUDE_PATH  ${GSL_INCLUDE_PATH})

#export vars
set (THIRD_PARTY_INCLUDE_PATH  ${THIRD_PARTY_INCLUDE_PATH} PARENT_SCOPE)
set (TEST_THIRD_PARTY_INCLUDE_PATH  ${TEST_THIRD_PARTY_INCLUDE_PATH} PARENT_SCOPE)
