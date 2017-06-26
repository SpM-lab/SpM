//=============================================================================
/*! _dsymatrix*double operator */
inline _dsymatrix operator*(const _dsymatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      mat.darray[j][i] *=d;
    }
  }
  
  return mat;
}

//=============================================================================
/*! dsymatrix/double operator */
inline _dsymatrix operator/(const _dsymatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      mat.darray[j][i] /=d;
    }
  }
  
  return mat;
}
