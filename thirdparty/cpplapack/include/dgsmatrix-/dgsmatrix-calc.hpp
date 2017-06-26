//=============================================================================
/*! return transposed dgsmatrix */
inline _dgsmatrix t(const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dgsmatrix newmat =mat;
  
  std::swap(newmat.m,newmat.n);
  std::swap(newmat.rows,newmat.cols);
  const std::vector<dcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<dcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    std::swap(it->i,it->j);
  }
  
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! search the index of element having the largest absolute value in 0-based numbering system */
inline void idamax(CPPL_INT& i, CPPL_INT& j, const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  //////// exception ////////
  if(mat.data.size()==0){
    WARNING_REPORT;
    std::cerr << "The dgsmatrix is a zero matrix." << std::endl;
    return;
  }
  
  //////// find ////////
  std::vector<dcomponent>::const_iterator itx(mat.data.begin());
  double vmax =0.;
  
  const std::vector<dcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    if( vmax < fabs(it->v) ){
      vmax =fabs(it->v);
      itx =it;
    }
  }
  
  i =itx->i;
  j =itx->j;
}

//=============================================================================
/*! return its largest absolute value */
inline double damax(const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  //////// exception ////////
  if(mat.data.size()==0){
    return 0.;
  }
  
  //////// find ////////
  std::vector<dcomponent>::const_iterator itx(mat.data.begin());
  double vmax =0.;
  
  const std::vector<dcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    if( vmax < fabs(it->v) ){
      vmax =fabs(it->v);
      itx =it;
    }
  }

  return itx->v;
}
