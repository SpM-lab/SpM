//=============================================================================
/*! double*_dgsmatrix operator */
inline _dgsmatrix operator*(const double& d, const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v *= d;
  }
  
  return mat;
}
