# find_package(Threads REQUIRED)
find_package(benchmark REQUIRED)

#~~~~~~~~~~~~~~~~

# message(STATUS "Test: ${LIBRARY_INCLUDE_PATH} ")
include_directories (${LIBRARY_INCLUDE_PATH} ${TEST_THIRD_PARTY_INCLUDE_PATH} ${xtensor_INCLUDE_DIRS})

file(GLOB_RECURSE ALL_BENCH_CPP *.cpp)

foreach(ONE_BENCH_CPP ${ALL_BENCH_CPP})

   get_filename_component(ONE_BENCH_EXEC ${ONE_BENCH_CPP} NAME_WE)

   # Avoid name collision 
   set(TARGET_NAME Bench_${ONE_BENCH_EXEC})

   add_executable(${TARGET_NAME} ${ONE_BENCH_CPP})
   set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME ${ONE_BENCH_EXEC}) 
   target_link_libraries(${TARGET_NAME} 

   # If required, you can add your project library here
# ...
    ${LIB_NAME} ${LIBS} fmt::fmt Threads::Threads benchmark::benchmark
   )

   # If you want to run benchmarks with the "make test" command, uncomment me
   add_test(${TARGET_NAME} ${ONE_BENCH_EXEC})
endforeach()
