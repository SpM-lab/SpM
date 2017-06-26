//=============================================================================
/*! +_dgsmatrix operator */
inline const _dgsmatrix& operator+(const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -_dgsmatrix operator */
inline _dgsmatrix operator-(const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v = -it->v;
  }
  
  return mat;
}
