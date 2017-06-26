//=============================================================================
/*! cast to _zgsmatrix */
inline _zgsmatrix _dgsmatrix::to_zgsmatrix() const
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat(m,n,CPPL_INT(data.size()));
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat.put(it->i, it->j, comple(it->v,0.));
  }
  
  destroy();
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to _dgematrix */
inline _dgematrix _dgsmatrix::to_dgematrix() const
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(m,n);
  newmat.zero();
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat(it->i,it->j) = it->v;
  }
  
  destroy();
  return _(newmat);
}
