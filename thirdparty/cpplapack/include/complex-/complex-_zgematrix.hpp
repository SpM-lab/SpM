//=============================================================================
/*! comple*_zgematrix operator */
inline _zgematrix operator*(const comple& d, const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  zscal_(&size, &d, mat.array, &inc);
  return mat;
}
