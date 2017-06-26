//=============================================================================
/*! _dgematrix*double operator */
inline _dgematrix operator*(const _dgematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT mn =mat.m*mat.n;
  CPPL_INT inc =1;
  dscal_(&mn, &d, mat.array, &inc);
  return mat;
}

//=============================================================================
/*! _dgematrix/double operator */
inline _dgematrix operator/(const _dgematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT mn =mat.m*mat.n;
  double dinv =1./d;
  CPPL_INT inc =1;
  dscal_(&mn, &dinv, mat.array, &inc);
  return mat;
}
