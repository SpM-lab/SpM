//=============================================================================
/*! double*_dgematrix operator */
inline _dgematrix operator*(const double& d, const _dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  dscal_(&size, &d, mat.array, &inc);
  return mat;
}
