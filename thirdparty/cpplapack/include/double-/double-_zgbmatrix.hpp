//=============================================================================
/*! double*_zgbmatrix operator */
inline _zgbmatrix operator*(const double& d, const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, mat.array, &inc);
  return mat;
}
