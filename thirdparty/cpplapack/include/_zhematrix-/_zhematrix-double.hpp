//=============================================================================
/*! _zhematrix*double operator */
inline _zhematrix operator*(const _zhematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.n*mat.n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, mat.array, &inc);
  return mat;
}

//=============================================================================
/*! _zhematrix/double operator */
inline _zhematrix operator/(const _zhematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.n*mat.n;
  double dinv =1./d;
  CPPL_INT inc =1;
  zdscal_(&size, &dinv, mat.array, &inc);
  return mat;
}
