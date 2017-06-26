//=============================================================================
/*! cast to _zgematrix */
inline _zgematrix _dgematrix::to_zgematrix() const
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(m,n);
  for(CPPL_INT i=0; i<m*n; i++){
    newmat.array[i] =comple(array[i],0.0);
  }
  
  destroy();
  return _(newmat);
}
