//=============================================================================
/*! comple*zgbmatrix operator */
inline _zgbmatrix operator*(const comple& d, const zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(newmat.kl+newmat.ku+1)*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =d*mat.array[i];
  }
  
  return _(newmat);
}
