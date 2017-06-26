//=============================================================================
/*! _drovector*double operator */
inline _drovector operator*(const _drovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  dscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}

//=============================================================================
/*! _drovector/double operator */
inline _drovector operator/(const _drovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  double dinv =1./d;
  CPPL_INT inc =1;
  dscal_(&vec.l, &dinv, vec.array, &inc);
  return vec;
}
