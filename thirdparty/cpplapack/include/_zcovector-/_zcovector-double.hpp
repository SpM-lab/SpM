//=============================================================================
/*! _zcovector*double operator */
inline _zcovector operator*(const _zcovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zdscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}

//=============================================================================
/*! _zcovector/double operator */
inline _zcovector operator/(const _zcovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  double dinv =1./d;
  CPPL_INT inc =1;
  zdscal_(&vec.l, &dinv, vec.array, &inc);
  return vec;
}
