//=============================================================================
/*! double*_dgbmatrix operator */
inline _dgbmatrix operator*(const double& d, const _dgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  CPPL_INT inc =1;
  dscal_(&size, &d, mat.array, &inc);
  return mat;
}
