//=============================================================================
/*! convert to _zgematrix */
inline _zgematrix zhsmatrix::to_zgematrix() const
{CPPL_VERBOSE_REPORT;
  zgematrix newmat( zgematrix(m,n).zero() );
  
  const std::vector<zcomponent>::const_iterator data_end =data.end();
  for(std::vector<zcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat(it->i, it->j) =it->v;
    newmat(it->j, it->i) =std::conj(it->v);
  }
  
  return _(newmat);
}

//=============================================================================
/*! convert to _zhematrix */
inline _zhematrix zhsmatrix::to_zhematrix() const
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(n);
  newmat.zero();
  
  const std::vector<zcomponent>::const_iterator data_end =data.end();
  for(std::vector<zcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat(it->i, it->j) =it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! convert to _zgsmatrix */
inline _zgsmatrix zhsmatrix::to_zgsmatrix() const
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat(m,n);
  newmat.zero();
  
  const std::vector<zcomponent>::const_iterator data_end =data.end();
  for(std::vector<zcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat.put(it->i, it->j, it->v);
    if(it->i!=it->j){
      newmat.put(it->j, it->i, std::conj(it->v));
    }
  }
  
  return _(newmat);
}
