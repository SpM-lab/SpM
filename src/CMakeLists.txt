#include directories
#include_directories(../include)

add_library(spm common.cpp SVD_matrix.cpp admm_svd.cpp spm_core.cpp kernel.cpp set_initial.cpp)
#execable
add_executable(SpM.out G2spectrum.cpp )
target_link_libraries(SpM.out spm ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
