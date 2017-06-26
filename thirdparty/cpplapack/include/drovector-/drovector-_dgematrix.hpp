//=============================================================================
/*! drovector*_dgematrix operator */
inline _drovector operator*(const drovector& vec, const _dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vec.l!=mat.m){
    ERROR_REPORT;
    std::cerr << "These vector and matrix can not make a product." << std::endl
              << "Your input was (" << vec.l << ") * (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  drovector newvec(mat.n);
  char trans ='T';
  double alpha =1.;
  CPPL_INT inc =1;
  double beta =0.;
  
  dgemv_( &trans, &mat.m, &mat.n, &alpha, mat.array, &mat.m, vec.array, &inc, &beta, newvec.array, &inc );
  
  mat.destroy();
  return _(newvec);
}
