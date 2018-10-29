#include directories
#include_directories(../include)

# SPM library
add_library(spm common.cpp SVD_matrix.cpp admm_svd.cpp spm_core.cpp kernel.cpp set_initial.cpp)

# GF library
add_library(gf Gf.cpp fft.cpp pade.cpp)
target_link_libraries(gf ${FFTW_LIBRARIES})

#execable
add_executable(SpM.out G2spectrum.cpp )
target_link_libraries(SpM.out spm gf ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
