cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

project(sparse CXX C)

# add common modules from ../common/cmake
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -D_DEBUG -g -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wconversion -Wdisabled-optimization -Wendif-labels -Wfloat-equal -Winit-self -Winline -Wlogical-op -Wmissing-include-dirs -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wpacked -Wpointer-arith -Wredundant-decls -Wshadow -Wsign-promo -Wswitch-default -Wswitch-enum -Wunsafe-loop-optimizations -Wvariadic-macros -Wwrite-strings ")

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3")

#include directories
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/cpplapack/include)

#find_package(OpenMP)
#if(OPENMP_FOUND)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif(OPENMP_FOUND)

if(USE_MKL)
  if(NOT MKLROOT)
    set(MKLROOT $ENV{MKLROOT})
  endif()
  if(NOT MKLROOT)
    message(FATAL "MKLROOT is not set")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
  include_directories(${MKLROOT}/include/fftw)
else()
  #Find BLAS
  find_package(BLAS REQUIRED)

  #Find LAPACK
  find_package(LAPACK REQUIRED)

  #Find FFTW3
  find_package(FFTW REQUIRED)
  include_directories(${FFTW_INCLUDE_DIRS})
endif()

# Build and enable tests
# testing setup
# enable_testing() must be called in the top-level CMakeLists.txt before any add_subdirectory() is called.
option(Testing "Enable testing" ON)

if (Testing)
  enable_testing(test)
  add_subdirectory(test/c++)
  add_subdirectory(test/python)
endif()

#Find SWIG
#find_package(SWIG REQUIRED)
#include(${SWIG_USE_FILE})

# search for python
#find_package(PythonLibs REQUIRED)
# NOTE: we need Python, but not necessarily the same version as the libraries (as of now!)
#find_package(PythonInterp REQUIRED)

# check for version mismatch.
#if(PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND AND NOT(PYTHON_VERSION_STRING VERSION_EQUAL PYTHONLIBS_VERSION_STRING))
#message(WARNING "Python executable (\"${PYTHON_VERSION_STRING}\") and library (\"${PYTHONLIBS_VERSION_STRING}\") version mismatch.")
#endif()

#include_directories(${PYTHON_INCLUDE_DIRS})

# find numpy
#find_package(NumPy REQUIRED)
#include_directories(${NUMPY_INCLUDE_DIRS})

# CPPLAPACK
#include_directories(${CPPLAPACK_INCLUDE_DIR})

add_subdirectory(c++)
