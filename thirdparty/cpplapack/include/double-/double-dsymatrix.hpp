//=============================================================================
/*! double*dsymatrix operator */
inline _dsymatrix operator*(const double& d, const dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  dsymatrix newmat(mat.n);
  
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      newmat.darray[j][i] =d*mat.darray[j][i];
    }
  }
  
  return _(newmat);
}
