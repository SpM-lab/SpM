//=============================================================================
/*! zcovector*=comple operator */
inline zcovector& zcovector::operator*=(const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zscal_(&l, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! zcovector/=comple operator */
inline zcovector& zcovector::operator/=(const comple& d)
{CPPL_VERBOSE_REPORT;
  comple dinv =1./d;
  CPPL_INT inc =1;
  zscal_(&l, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zcovector*comple operator */
inline _zcovector operator*(const zcovector& vec, const comple& d)
{CPPL_VERBOSE_REPORT;
  zcovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]*d;
  }
  
  return _(newvec);
}

//=============================================================================
/*! zcovector/comple operator */
inline _zcovector operator/(const zcovector& vec, const comple& d)
{CPPL_VERBOSE_REPORT;
  zcovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]/d;
  }
  
  return _(newvec);
}
