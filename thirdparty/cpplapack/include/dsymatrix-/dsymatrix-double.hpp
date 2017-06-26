//=============================================================================
/*! dsymatrix*=double operator */
inline dsymatrix& dsymatrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      darray[j][i] *=d;
    }
  }
  
  return *this;
}

//=============================================================================
/*! dsymatrix/=double operator */
inline dsymatrix& dsymatrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      darray[j][i] /=d;
    }
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dsymatrix*double operator */
inline _dsymatrix operator*(const dsymatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dsymatrix newmat(mat.n);
  
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      newmat.darray[j][i] =mat.darray[j][i]*d;
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! dsymatrix/double operator */
inline _dsymatrix operator/(const dsymatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dsymatrix newmat(mat.n);
  
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      newmat.darray[j][i] =mat.darray[j][i]/d;
    }
  }
  
  return _(newmat);
}
