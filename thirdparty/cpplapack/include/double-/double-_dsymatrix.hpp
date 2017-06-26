//=============================================================================
/*! double*_dsymatrix operator */
inline _dsymatrix operator*(const double& d, const _dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      mat.darray[j][i] *=d;
    }
  }
  
  return mat;
}
