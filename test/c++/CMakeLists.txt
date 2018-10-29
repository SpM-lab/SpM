cmake_minimum_required(VERSION 2.8.12)

function(add_gtest test)
    set(source "${ARGV1}/${test}.cpp")
    #set(gtest_src "${ARGV1}/gtest_main.cc;${ARGV1}/gtest-all.cc")
    
    add_executable(${test} ${source} ${header_files})
    target_link_libraries(${test} GTest::GTest GTest::Main)
    target_include_directories(${test} PRIVATE ${HDF5_INCLUDE_DIR})
    add_test(NAME ${test} COMMAND ${test} ${test_xml_output})
endfunction(add_gtest)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(../../c++/include)
include_directories(./)

# test for SPM library


# test for GF library
set(unittest_src_gf
test_fft
test_pade
)

foreach(test_c++ ${unittest_src_gf})
  add_gtest(${test_c++} ".")
  target_link_libraries(${test_c++} gf)
endforeach(test_c++)

add_subdirectory(data)
