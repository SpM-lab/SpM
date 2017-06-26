//=============================================================================
/*! cast to _zgbmatrix */
inline _zgbmatrix _dgbmatrix::to_zgbmatrix() const
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(m,n,kl,ku);
  
  for(CPPL_INT i=0; i<(kl+ku+1)*n; i++){
    newmat.array[i] =comple(array[i],0.0);
  }
  
  destroy();
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to _dgematrix */
inline _dgematrix _dgbmatrix::to_dgematrix() const
{CPPL_VERBOSE_REPORT;
  dgematrix newmat( dgematrix(m,n).zero() );
  
  for(CPPL_INT i=0; i<m; i++){
    const CPPL_INT jmax =std::min(n,i+ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-kl); j<jmax; j++){
      newmat(i,j) =(*this)(i,j);
    }
  }
  
  destroy();
  return _(newmat);
}
