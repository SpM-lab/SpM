//=============================================================================
/*! dgbmatrix*_dcovector operator */
inline _dcovector operator*(const dgbmatrix& mat, const _dcovector& vec)
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
  CPPL_INT lda =mat.kl+mat.ku+1;
  CPPL_INT inc =1;
  double beta =0.;

  dgbmv_( &trans, &mat.m, &mat.n, &mat.kl, &mat.ku, &alpha, mat.array, &lda, vec.array, &inc, &beta, newvec.array, &inc );
  
  vec.destroy();
  return _(newvec);
}
