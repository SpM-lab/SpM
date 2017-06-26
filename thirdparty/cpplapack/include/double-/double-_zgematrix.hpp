//=============================================================================
/*! double*_zgematrix operator */
inline _zgematrix operator*(const double& d, const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, mat.array, &inc);
  return mat;
}
