//=============================================================================
/*! _zrovector*zhematrix operator */
inline _zrovector operator*(const _zrovector& vec, const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vec.l!=mat.n){
    ERROR_REPORT;
    std::cerr << "These vector and matrix can not make a product." << std::endl
              << "Your input was (" << vec.l << ") * (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zrovector newvec(mat.n);
  char uplo ='l';
  comple alpha =comple(1.,0.);
  CPPL_INT inc =1;
  comple beta =comple(0.,0.);
  
  zhemv_( &uplo, &mat.n, &alpha, mat.array, &mat.n, vec.array, &inc, &beta, newvec.array, &inc );
  
  vec.destroy();
  return _(newvec);
}
