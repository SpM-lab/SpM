//=============================================================================
/*! zcovector*zrovector operator */
inline _zgematrix operator*(const zcovector& covec, const zrovector& rovec)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(covec.l, rovec.l);
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =covec(i)*rovec(j);
    }
  }
  
  return _(newmat);
}
