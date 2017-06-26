//=============================================================================
/*! convert to _zgematrix */
inline _zgematrix _zgsmatrix::to_zgematrix() const
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(m,n);
  newmat.zero();
  
  const std::vector<zcomponent>::const_iterator data_end =data.end();
  for(std::vector<zcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat(it->i,it->j) = it->v;
  }
  
  destroy();
  return _(newmat);
}
