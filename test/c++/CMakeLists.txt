cmake_minimum_required(VERSION 2.8.12)
include(EnableGtests) #defined in ./cmake
include_directories(../../c++/include)

# test for SPM library


# test for GF library
set(unittest_src_gf
test_fft_fermion
test_fft_boson
test_pade_fermion
test_pade_boson
)

foreach(test_c++ ${unittest_src_gf})
  add_gtest(${test_c++})
  target_link_libraries(${test_c++} gf)
endforeach(test_c++)

add_subdirectory(data)
