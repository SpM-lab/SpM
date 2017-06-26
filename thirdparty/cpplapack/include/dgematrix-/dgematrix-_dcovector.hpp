//=============================================================================
/*! dgematrix*_dcovector operator */
inline _dcovector operator*(const dgematrix& mat, const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.n!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector can not make a product." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ") * (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dcovector newvec(mat.m);
  char trans ='n';
  double alpha =1.;
  CPPL_INT inc =1;
  double beta =0.;
  
  dgemv_( &trans, &mat.m, &mat.n, &alpha, mat.array, &mat.m, vec.array, &inc, &beta, newvec.array, &inc );
  
  vec.destroy();
  return _(newvec);
}
