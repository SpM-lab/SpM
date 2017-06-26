//=============================================================================
/*! _dsymatrix*_dcovector operator */
inline _dcovector operator*(const _dsymatrix& mat, const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.n!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector can not make a product." << std::endl
              << "Your input was (" << mat.n << "x" << mat.n << ") * (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dcovector newvec(mat.n);
  char uplo ='l';
  double alpha =1.;
  CPPL_INT inc =1;
  double beta =0.;
  
  dsymv_( &uplo, &mat.n, &alpha, mat.array, &mat.n, vec.array, &inc, &beta, newvec.array, &inc );
  
  mat.destroy();
  vec.destroy();
  return _(newvec);
}
