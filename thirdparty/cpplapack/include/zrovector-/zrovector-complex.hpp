//=============================================================================
/*! zrovector*=comple operator */
inline zrovector& zrovector::operator*=(const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zscal_(&l, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! zrovector/=comple operator */
inline zrovector& zrovector::operator/=(const comple& d)
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
/*! zrovector*comple operator */
inline _zrovector operator*(const zrovector& vec, const comple& d)
{CPPL_VERBOSE_REPORT;
  zrovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]*d;
  }
  
  return _(newvec);
}

//=============================================================================
/*! zrovector/comple operator */
inline _zrovector operator/(const zrovector& vec, const comple& d)
{CPPL_VERBOSE_REPORT;
  zrovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]/d;
  }
  
  return _(newvec);
}
