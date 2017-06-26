//=============================================================================
/*! zgematrix*=comple operator */
inline zgematrix& zgematrix::operator*=(const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =m*n;
  CPPL_INT inc =1;
  zscal_(&size, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! zgematrix/=comple operator */
inline zgematrix& zgematrix::operator/=(const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =m*n;
  comple dinv =1./d;
  CPPL_INT inc =1;
  zscal_(&size, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix*comple operator */
inline _zgematrix operator*(const zgematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.m, mat.n);
  
  const CPPL_INT size =mat.m*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]*d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zgematrix/comple operator */
inline _zgematrix operator/(const zgematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.m, mat.n);
  
  const CPPL_INT size =mat.m*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]/d;
  }
  
  return _(newmat);
}
