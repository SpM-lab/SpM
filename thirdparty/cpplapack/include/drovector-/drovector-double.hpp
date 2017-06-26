//=============================================================================
/*! drovector*=double operator */
inline drovector& drovector::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  dscal_(&l, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! drovector/=double operator */
inline drovector& drovector::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  double dinv =1./d;
  CPPL_INT inc =1;
  dscal_(&l, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! drovector*double operator */
inline _drovector operator*(const drovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  drovector newvec(vec.l);
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]*d;
  }
  
  return _(newvec);
}

//=============================================================================
/*! drovector/double operator */
inline _drovector operator/(const drovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  drovector newvec(vec.l);
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =vec.array[i]/d;
  }
  
  return _(newvec);
}
