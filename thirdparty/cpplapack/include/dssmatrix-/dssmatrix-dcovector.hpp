//=============================================================================
/*! dssmatrix*dcovector operator */
inline _dcovector operator*(const dssmatrix& mat, const dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.n!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector can not make a product." << std::endl
              << "Your input was (" << mat.n << "x" << mat.n << ") * (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  dcovector newvec(mat.n);
  newvec.zero();
  
  const std::vector<dcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    newvec(it->i) +=it->v*vec(it->j);
    if(it->i!=it->j){
      newvec(it->j) +=it->v*vec(it->i);
    }
  }
  
  return _(newvec);
}
