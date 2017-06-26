//=============================================================================
/*! convert to _zgematrix */
inline _zgematrix _zhematrix::to_zgematrix() const
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(n,n);
  
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<n; j++){
      newmat(i,j) =(*this)(i,j);
    }
  }
  
  destroy();
  return _(newmat);
}
