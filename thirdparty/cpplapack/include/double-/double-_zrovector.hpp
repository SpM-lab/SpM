//=============================================================================
/*! double*_zrovector operator */
inline _zrovector operator*(const double& d, const _zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zdscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}
