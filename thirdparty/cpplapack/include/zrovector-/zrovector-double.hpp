//=============================================================================
/*! zrovector*=double operator */
inline zrovector& zrovector::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  zdscal_(&l, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! zrovector/=double operator */
inline zrovector& zrovector::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  double dinv =1./d;
  CPPL_INT inc =1;
  zdscal_(&l, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector*double operator */
inline _zrovector operator*(const zrovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  zrovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]*d;
  }
  
  return _(newvec);
}

//=============================================================================
/*! zrovector/double operator */
inline _zrovector operator/(const zrovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  zrovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]/d;
  }
  
  return _(newvec);
}
