//=============================================================================
/*! comple*_zrovector operator */
inline _zrovector operator*(const comple& d, const _zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}
