//=============================================================================
/*! double*dgbmatrix operator */
inline _dgbmatrix operator*(const double& d, const dgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =d*mat.array[i];
  }
  
  return _(newmat);
}
