//=============================================================================
/*! double*dgematrix operator */
inline _dgematrix operator*(const double& d, const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(mat.m, mat.n);
  
  const CPPL_INT size =mat.m*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =d*mat.array[i];
  }
  
  return _(newmat);
}
