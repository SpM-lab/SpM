//=============================================================================
/*! complex*zgsmatrix operator */
inline _zgsmatrix operator*(const comple& d, const zgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat =mat;
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v *= d;
  }
  
  return _(newmat);
}
