//=============================================================================
/*! _zrovector*comple operator */
inline _zrovector operator*(const _zrovector& vec, const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zscal_(&vec.l, &d, vec.array, &inc);
  return vec;
}

//=============================================================================
/*! _zrovector/comple operator */
inline _zrovector operator/(const _zrovector& vec, const comple& d)
{CPPL_VERBOSE_REPORT;
  comple dinv =1./d;
  CPPL_INT inc =1;
  zscal_(&vec.l, &dinv, vec.array, &inc);
  return vec;
}
