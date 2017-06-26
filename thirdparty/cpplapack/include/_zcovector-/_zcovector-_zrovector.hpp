//=============================================================================
/*! _zcovector*_zrovector operator */
inline _zgematrix operator*(const _zcovector& covec, const _zrovector& rovec)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(covec.l, rovec.l);
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =covec(i)*rovec(j);
    }
  }
  
  covec.destroy();
  rovec.destroy();
  return _(newmat);
}
