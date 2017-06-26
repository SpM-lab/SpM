//=============================================================================
/*! _dssmatrix*double operator */
inline _dssmatrix operator*(const _dssmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v *=d;
  }
  
  return mat;
}

//=============================================================================
/*! _dssmatrix/double operator */
inline _dssmatrix operator/(const _dssmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v /=d;
  }
  
  return mat;
}
