//=============================================================================
/*! comple*_zgbmatrix operator */
inline _zgbmatrix operator*(const comple& d, const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  CPPL_INT inc =1;
  zscal_(&size, &d, mat.array, &inc);
  return mat;
}
