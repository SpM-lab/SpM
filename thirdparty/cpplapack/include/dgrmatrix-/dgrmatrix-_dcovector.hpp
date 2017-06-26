//=============================================================================
/*! dgrmatrix*_dcovector operator */
inline _dcovector operator*(const dgrmatrix& mat, const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.n!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector can not make a product." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ") * (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dcovector VEC =vec;
  dcovector newvec =mat*VEC;
  return _(newvec);
}
