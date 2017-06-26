//=============================================================================
/*! _drovector*_dgbmatrix operator */
inline _drovector operator*(const _drovector& vec, const _dgbmatrix& mat)
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
  CPPL_INT lda =mat.kl+mat.ku+1;
  CPPL_INT inc =1;
  double beta =0.;

  dgbmv_( &trans, &mat.m, &mat.n, &mat.kl, &mat.ku, &alpha, mat.array, &lda, vec.array, &inc, &beta, newvec.array, &inc );
  
  vec.destroy();
  mat.destroy();
  return _(newvec);
}
