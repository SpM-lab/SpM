//=============================================================================
/*! +_zgsmatrix operator */
inline const _zgsmatrix& operator+(const _zgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -_zgsmatrix operator */
inline _zgsmatrix operator-(const _zgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v = -it->v;
  }
  
  return mat;
}
