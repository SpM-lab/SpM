//=============================================================================
/*! comple*zgematrix operator */
inline _zgematrix operator*(const comple& d, const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.m, mat.n);
  
  const CPPL_INT size =mat.m*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =d*mat.array[i];
  }
  
  return _(newmat);
}
