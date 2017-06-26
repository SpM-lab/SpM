//=============================================================================
/*! cast to _zhematrix */
inline _zhematrix _dsymatrix::to_zhematrix() const
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(n);
  
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      newmat(i,j) =comple((*this)(i,j),0.0);
    }
  }
  
  destroy();
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to _dgematrix */
inline _dgematrix _dsymatrix::to_dgematrix() const
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(n,n);
  
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<n; j++){
      newmat(i,j) =(*this)(i,j);
    }
  }
  
  destroy();
  return _(newmat);
}
