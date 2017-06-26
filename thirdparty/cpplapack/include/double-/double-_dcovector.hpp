//=============================================================================
/*! double*_dcovector operator */
inline _dcovector operator*(const double& d, const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  dscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}
