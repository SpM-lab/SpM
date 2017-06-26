//=============================================================================
/*! cast to _zhsmatrix */
inline _zhsmatrix _dssmatrix::to_zhsmatrix() const
{CPPL_VERBOSE_REPORT;
  zhsmatrix newmat(n,CPPL_INT(data.size()));
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat.put(it->i, it->j, comple(it->v,0.0));
  }
  
  destroy();
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to _dgematrix */
inline _dgematrix _dssmatrix::to_dgematrix() const
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(m,n);
  newmat.zero();
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat(it->i, it->j) =it->v;
    newmat(it->j, it->i) =it->v;
  }
  
  destroy();
  return _(newmat);
}

//=============================================================================
/*! convert to _dsymatrix */
inline _dsymatrix _dssmatrix::to_dsymatrix() const
{CPPL_VERBOSE_REPORT;
  dsymatrix newmat(n);
  newmat.zero();
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat(it->i, it->j) =it->v;
  }
  
  destroy();
  return _(newmat);
}

//=============================================================================
/*! convert to _dgsmatrix */
inline _dgsmatrix _dssmatrix::to_dgsmatrix() const
{CPPL_VERBOSE_REPORT;
  dgsmatrix newmat(m,n,CPPL_INT(data.size()*2));
  newmat.zero();
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat.put(it->i, it->j, it->v);
    if(it->i!=it->j){
      newmat.put(it->j, it->i, it->v);
    }
  }
  
  destroy();
  return _(newmat);
}
