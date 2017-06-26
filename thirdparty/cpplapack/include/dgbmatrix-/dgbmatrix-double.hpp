//=============================================================================
/*! dgbmatrix*=double operator */
inline dgbmatrix& dgbmatrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(kl+ku+1)*n;
  CPPL_INT inc =1;
  dscal_(&size, &d, array, &inc);
  return *this;
}

//=============================================================================
/*! dgbmatrix/=double operator */
inline dgbmatrix& dgbmatrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(kl+ku+1)*n;
  double dinv =1./d;
  CPPL_INT inc =1;
  dscal_(&size, &dinv, array, &inc);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgbmatrix*double operator */
inline _dgbmatrix operator*(const dgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(newmat.kl+newmat.ku+1)*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]*d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgbmatrix/double operator */
inline _dgbmatrix operator/(const dgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgbmatrix newmat(mat.m, mat.n, mat.kl, mat.ku);
  
  const CPPL_INT size =(newmat.kl+newmat.ku+1)*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]/d;
  }
  
  return _(newmat);
}
