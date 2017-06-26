//=============================================================================
/*! _zgematrix*double operator */
inline _zgematrix operator*(const _zgematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, mat.array, &inc);
  return mat;
}

//=============================================================================
/*! _zgematrix/double operator */
inline _zgematrix operator/(const _zgematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  double dinv =1./d;
  CPPL_INT inc =1;
  zdscal_(&size, &dinv, mat.array, &inc);
  return mat;
}
