//=============================================================================
/*! double*_zhsmatrix operator */
inline _zhsmatrix operator*(const double& d, const _zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::iterator it=mat.data.begin(); it!=mat_data_end; it++){
    it->v *=d;
  }
  
  return mat;
}
