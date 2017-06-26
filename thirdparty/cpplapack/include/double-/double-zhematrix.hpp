//=============================================================================
/*! double*zhematrix operator */
inline _zhematrix operator*(const double& d, const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(mat.n);
  
  const CPPL_INT size =mat.n*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =d*mat.array[i];
  }
  
  return _(newmat);
}
