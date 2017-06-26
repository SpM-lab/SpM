//=============================================================================
/*! drovector*_dsymatrix operator */
inline _drovector operator*(const drovector& vec, const _dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vec.l!=mat.n){
    ERROR_REPORT;
    std::cerr << "These vector and matrix can not make a product." << std::endl
              << "Your input was (" << vec.l << ") * (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  drovector newvec(mat.n);
  char uplo ='l';
  double alpha =1.;
  CPPL_INT inc =1;
  double beta =0.;
  
  dsymv_( &uplo, &mat.n, &alpha, mat.array, &mat.n, vec.array, &inc, &beta, newvec.array, &inc );
  
  mat.destroy();
  return _(newvec);
}
