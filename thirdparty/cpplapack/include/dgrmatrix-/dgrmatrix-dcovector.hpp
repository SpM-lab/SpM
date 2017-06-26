//=============================================================================
/*! dgrmatrix*dcovector operator */
inline _dcovector operator*(const dgrmatrix& mat, const dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.n!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector can not make a product." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ") * (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
#ifdef  __INTEL_COMPILER
  dcovector newvec(mat.m);
  char transa ='N';
  MKL_INT m =MKL_INT(mat.m);
  double* a =const_cast<double*>(&mat.a[0]);
  MKL_INT* ia =const_cast<MKL_INT*>(&mat.ia[0]);
  MKL_INT* ja =const_cast<MKL_INT*>(&mat.ja[0]);
  MKL_DCSRGEMV(&transa, &m, a, ia, ja, vec.array, newvec.array);
  return _(newvec);
  
  
#else //__INTEL_COMPILER is not defined
  dcovector newvec(mat.m);  
#pragma omp parallel for
  for(CPPL_INT i=0; i<mat.m; i++){
    double sum =0.;
    int k_beg =mat.ia[i]-1;
    int k_end =mat.ia[i+1]-1;
    for(int k=k_beg; k<k_end; k++){
      int j =mat.ja[k]-1;
      sum += mat.a[k] * vec(j);
    }
    newvec(i) =sum;
  }
  return _(newvec);
#endif//__INTEL_COMPILER
}
