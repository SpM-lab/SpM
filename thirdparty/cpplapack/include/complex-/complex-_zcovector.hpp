//=============================================================================
/*! comple*_zcovector operator */
inline _zcovector operator*(const comple& d, const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zscal_(&vec.l, &d, vec.array, &inc);  
  return vec;
}
