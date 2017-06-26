//=============================================================================
/*! zgbmatrix*=double operator */
inline zgbmatrix& zgbmatrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(kl+ku+1)*n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! zgbmatrix/=double operator */
inline zgbmatrix& zgbmatrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(kl+ku+1)*n;
  double dinv =1./d;
  CPPL_INT inc =1;
  zdscal_(&size, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgbmatrix*double operator */
inline _zgbmatrix operator*(const zgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(newmat.kl+newmat.ku+1)*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]*d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zgbmatrix/double operator */
inline _zgbmatrix operator/(const zgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(newmat.kl+newmat.ku+1)*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]/d;
  }
  
  return _(newmat);
}
