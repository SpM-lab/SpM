//=============================================================================
/*! _dgsmatrix+dgsmatrix operator */
inline _dgsmatrix operator+(const _dgsmatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgsmatrix newmat(matA);
  
  const std::vector<dcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<dcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    newmat(it->i,it->j) += it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! _dgsmatrix-dgsmatrix operator */
inline _dgsmatrix operator-(const _dgsmatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgsmatrix newmat(matA);
  
  const std::vector<dcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<dcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    newmat(it->i,it->j) -= it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! _dgsmatrix*dgsmatrix operator */
inline _dgsmatrix operator*(const _dgsmatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgsmatrix newmat( matA.m, matB.n );
  
  const std::vector<dcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<dcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    CPPL_INT k =it->j;
    const std::vector<CPPL_INT>::const_iterator matB_rows_k_end =matB.rows[k].end();
    for(std::vector<CPPL_INT>::const_iterator p=matB.rows[k].begin(); p!=matB_rows_k_end; p++){
      newmat(it->i,matB.data[*p].j) += it->v*matB.data[*p].v;
    }
  }
  
  matA.destroy();
  return _(newmat);
}
