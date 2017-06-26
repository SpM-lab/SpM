//=============================================================================
/*! zhematrix*=double operator */
inline zhematrix& zhematrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =n*n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! zhematrix/=double operator */
inline zhematrix& zhematrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =n*n;
  double dinv =1./d;
  CPPL_INT inc =1;
  zdscal_(&size, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix*double operator */
inline _zhematrix operator*(const zhematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(mat.n);
  
  const CPPL_INT size =mat.n*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]*d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zhematrix/double operator */
inline _zhematrix operator/(const zhematrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(mat.n);
  
  const CPPL_INT size =mat.n*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]/d;
  }
  
  return _(newmat);
}
