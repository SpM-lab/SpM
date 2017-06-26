//=============================================================================
/*! double*zhsmatrix operator */
inline _zhsmatrix operator*(const double& d, const zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  zhsmatrix newmat =mat;
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v *=d;
  }
  
  return _(newmat);
}
