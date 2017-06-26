//=============================================================================
/*! _dcovector*double operator */
inline _dcovector operator*(const _dcovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  dscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}

//=============================================================================
/*! _dcovector/double operator */
inline _dcovector operator/(const _dcovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  double dinv =1./d;
  dscal_(&vec.l, &dinv, vec.array, &inc);
  return vec;
}
