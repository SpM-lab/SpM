//=============================================================================
/*! convert to _zgematrix */
inline _zgematrix _zgbmatrix::to_zgematrix() const
{CPPL_VERBOSE_REPORT;
  zgematrix newmat( zgematrix(m,n).zero() );
  
  for(CPPL_INT i=0; i<m; i++){
    const CPPL_INT jmax =std::min(n,i+ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-kl); j<jmax; j++){
      newmat(i,j) =(*this)(i,j);
    }
  }
  
  destroy();
  return _(newmat);
}
