//=============================================================================
/*! zrovector*zgematrix operator */
inline _zrovector operator*(const zrovector& vec, const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vec.l!=mat.m){
    ERROR_REPORT;
    std::cerr << "These vector and matrix can not make a product." << std::endl
              << "Your input was (" << vec.l << ") * (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zrovector newvec(mat.n);
  char trans ='T';
  comple alpha =comple(1.,0.);
  CPPL_INT inc =1;
  comple beta =comple(0.,0.);
  
  zgemv_( &trans, &mat.m, &mat.n, &alpha, mat.array, &mat.m, vec.array, &inc, &beta, newvec.array, &inc );
  
  return _(newvec);
}
