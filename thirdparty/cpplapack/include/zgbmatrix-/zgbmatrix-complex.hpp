//=============================================================================
/*! zgbmatrix*=comple operator */
inline zgbmatrix& zgbmatrix::operator*=(const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(kl+ku+1)*n;
  CPPL_INT inc =1;
  zscal_(&size, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! zgbmatrix/=comple operator */
inline zgbmatrix& zgbmatrix::operator/=(const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(kl+ku+1)*n;
  comple dinv =1./d;
  CPPL_INT inc =1;
  zscal_(&size, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgbmatrix*comple operator */
inline _zgbmatrix operator*(const zgbmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(newmat.kl+newmat.ku+1)*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]*d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zgbmatrix/comple operator */
inline _zgbmatrix operator/(const zgbmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(newmat.kl+newmat.ku+1)*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]/d;
  }
  
  return _(newmat);
}
