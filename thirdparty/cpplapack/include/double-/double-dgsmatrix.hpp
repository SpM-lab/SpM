//=============================================================================
/*! double*dgsmatrix operator */
inline _dgsmatrix operator*(const double& d, const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dgsmatrix newmat =mat;
  
  const std::vector<dcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<dcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v *= d;
  }
  
  return _(newmat);
}
