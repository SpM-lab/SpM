//=============================================================================
/*! return transposed dssmatrix */
inline _dssmatrix t(const dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  WARNING_REPORT;
  std::cerr << "This function call has no effect since the matrix is symmetric." << std::endl;
#endif//CPPL_DEBUG

  dssmatrix newmat(mat);
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! search the index of element having the largest absolute value
  in 0-based numbering system */
inline void idamax(CPPL_INT& i, CPPL_INT& j, const dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
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
inline double damax(const dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
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
