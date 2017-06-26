//=============================================================================
/*! double*_zcovector operator */
inline _zcovector operator*(const double& d, const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zdscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}
