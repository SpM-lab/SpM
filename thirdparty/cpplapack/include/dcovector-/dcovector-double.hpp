//=============================================================================
/*! dcovector*=double operator */
inline dcovector& dcovector::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  dscal_(&l, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! dcovector/=double operator */
inline dcovector& dcovector::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  double dinv =1./d;
  dscal_(&l, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector*double operator */
inline _dcovector operator*(const dcovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  dcovector newvec(vec.l);
  for(CPPL_INT i=0; i<vec.l; i++){ newvec.array[i] =vec.array[i]*d; }
  
  return _(newvec);
}

//=============================================================================
/*! dcovector/double operator */
inline _dcovector operator/(const dcovector& vec, const double& d)
{CPPL_VERBOSE_REPORT;
  dcovector newvec(vec.l);
  for(CPPL_INT i=0; i<vec.l; i++){ newvec.array[i] =vec.array[i]/d; }
  
  return _(newvec);
}
