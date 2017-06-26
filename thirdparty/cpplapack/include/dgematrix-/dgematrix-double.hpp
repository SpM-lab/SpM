//=============================================================================
/*! dgematrix*=double operator */
inline dgematrix& dgematrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT mn =m*n;
  CPPL_INT inc =1;
  dscal_(&mn, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! dgematrix/=double operator */
inline dgematrix& dgematrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT mn =m*n;
  double dinv =1./d;
  CPPL_INT inc =1;
  dscal_(&mn, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgematrix*double operator */
inline _dgematrix operator*(const dgematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(mat.m, mat.n);
  
  const CPPL_INT mn =mat.m*mat.n;
  for(CPPL_INT i=0; i<mn; i++){
    newmat.array[i] =mat.array[i]*d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgematrix/double operator */
inline _dgematrix operator/(const dgematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(mat.m, mat.n);
  
  const CPPL_INT mn =mat.m*mat.n;
  for(CPPL_INT i=0; i<mn; i++){
    newmat.array[i] =mat.array[i]/d;
  }
  
  return _(newmat);
}
