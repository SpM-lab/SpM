//=============================================================================
/*! _zgsmatrix+zhematrix operator */
inline _zgematrix operator+(const _zgsmatrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.m!=matB.n || matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat( matB.to_zgematrix() );
  
  const std::vector<zcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<zcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    newmat(it->i,it->j) += it->v;
  }
  
  matA.destroy();
  return _(newmat);
}

//=============================================================================
/*! _zgsmatrix-zhematrix operator */
inline _zgematrix operator-(const _zgsmatrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.m!=matB.n || matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// shallow copy to zgematrix ////
  zgematrix newmat( (-matB).to_zgematrix() );
  
  //// add ////
  const std::vector<zcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<zcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    newmat(it->i,it->j) += it->v;
  }
  
  //// return ////
  matA.destroy();
  return _(newmat);
}

//=============================================================================
/*! _zgsmatrix*zhematrix operator */
inline _zgematrix operator*(const _zgsmatrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat(matA.m, matB.n);
  newmat.zero();
  
  const std::vector<zcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<zcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    for(CPPL_INT i=0; i<matB.n; i++){
      newmat(it->i,i) += it->v*matB(it->j,i);
    }
  }
  
  matA.destroy();
  return _(newmat);
}
