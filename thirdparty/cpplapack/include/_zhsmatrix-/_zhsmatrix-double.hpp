//=============================================================================
/*! _zhsmatrix*double operator */
inline _zhsmatrix operator*(const _zhsmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v *=d;
  }
  
  return mat;
}

//=============================================================================
/*! _zhsmatrix/double operator */
inline _zhsmatrix operator/(const _zhsmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v /=d;
  }
  
  return mat;
}
