//=============================================================================
/*! double*_drovector operator */
inline _drovector operator*(const double& d, const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  dscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}
