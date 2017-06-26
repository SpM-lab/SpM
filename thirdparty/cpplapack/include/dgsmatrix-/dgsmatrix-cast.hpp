//=============================================================================
/*! cast to _zgsmatrix */
inline _zgsmatrix dgsmatrix::to_zgsmatrix() const
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat(m,n,CPPL_INT(data.size()));
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat.put(it->i, it->j, comple(it->v,0.));
  }
  
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to _dgematrix */
inline _dgematrix dgsmatrix::to_dgematrix() const
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(m,n);
  newmat.zero();
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    newmat(it->i,it->j) = it->v;
  }
  
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to dgrmatrix */
inline dgrmatrix dgsmatrix::to_dgrmatrix() const
{CPPL_VERBOSE_REPORT;
  //////// resize ////////
  dgrmatrix newmat;
  newmat.m =m;
  newmat.n =n;
  newmat.a.resize(data.size());
  newmat.ia.resize(m+1);
  newmat.ja.resize(data.size());
  
  //////// copy ////////
  newmat.ia[0] =1;//one-based
  CPPL_INT k=0;
  for(CPPL_INT i=0; i<m; i++){
    //// make map ////
    const std::vector<CPPL_INT>::const_iterator rows_i_end =rows[i].end();
    std::map<CPPL_INT,CPPL_INT> jc;
    for(std::vector<CPPL_INT>::const_iterator rit=rows[i].begin(); rit!=rows_i_end; rit++){
      jc.insert( std::make_pair(data[*rit].j, *rit) );
    }
    //// assign ////
    const std::map<CPPL_INT,CPPL_INT>::const_iterator jc_end =jc.end();
    for(std::map<CPPL_INT,CPPL_INT>::const_iterator jcit=jc.begin(); jcit!=jc_end; jcit++){
      newmat.a[k] =data[(*jcit).second].v;
      newmat.ja[k] =CPPL_INT((*jcit).first)+1;//one-based
      k++;
    }
    newmat.ia[i+1] =k+1;//one-based
  }
  
  return newmat;
}
