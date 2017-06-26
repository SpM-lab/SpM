//=============================================================================
/*! comple*zhematrix operator */
inline _zgematrix operator*(const comple& d, const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  mat.complete();
  zgematrix newmat(mat.n, mat.n);
  
  const CPPL_INT size =mat.n*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =d*mat.array[i];
  }
  
  return _(newmat);
}
