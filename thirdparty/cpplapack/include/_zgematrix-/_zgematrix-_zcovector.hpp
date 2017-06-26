//=============================================================================
/*! _zgematrix*_zcovector operator */
inline _zcovector operator*(const _zgematrix& mat, const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.n!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector can not make a product." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ") * (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zcovector newvec(mat.m);
  char trans ='n';
  comple alpha =comple(1.,0.);
  CPPL_INT inc =1;
  comple beta =comple(0.,0.);
  
  zgemv_( &trans, &mat.m, &mat.n, &alpha, mat.array, &mat.m, vec.array, &inc, &beta, newvec.array, &inc );
  
  mat.destroy();
  vec.destroy();
  return _(newvec);
}
