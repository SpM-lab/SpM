//=============================================================================
/*! cast to _zhematrix */
inline _zhematrix dsymatrix::to_zhematrix() const
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(n);
  
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      newmat(i,j) =comple((*this)(i,j),0.0);
    }
  }
  
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to _dgematrix */
inline _dgematrix dsymatrix::to_dgematrix() const
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(n,n);
  
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<n; j++){
      newmat(i,j) =(*this)(i,j);
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! convert to _dssmatrix */
inline _dssmatrix dsymatrix::to_dssmatrix(const double eps) const
{CPPL_VERBOSE_REPORT;
  dssmatrix newmat(n);
  
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      if( fabs((*this)(i,j))>eps ){
        newmat(i,j) =(*this)(i,j);
      }
    }
  }
  
  return _(newmat);
}
