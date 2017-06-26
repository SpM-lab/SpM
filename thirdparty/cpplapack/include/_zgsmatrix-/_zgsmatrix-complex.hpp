//=============================================================================
/*! _zgsmatrix*comple operator */
inline _zgsmatrix operator*(const _zgsmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v *=d;
  }
  
  return mat;
}

//=============================================================================
/*! _zgsmatrix/comple operator */
inline _zgsmatrix operator/(const _zgsmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v /=d;
  }
  
  return mat;
}
