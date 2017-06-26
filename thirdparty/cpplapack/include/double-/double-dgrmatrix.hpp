//=============================================================================
/*! double*dgrmatrix operator */
inline dgrmatrix operator*(const double& d, const dgrmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dgrmatrix newmat =mat;
  
  const size_t newmat_a_size =newmat.a.size();
  for(size_t k=0; k<newmat_a_size; k++){
    newmat.a[k] *= d;
  }
  
  return newmat;
}
