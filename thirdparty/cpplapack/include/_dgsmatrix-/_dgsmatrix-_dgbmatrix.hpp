//=============================================================================
/*! _dgsmatrix+_dgbmatrix operator */
inline _dgematrix operator+(const _dgsmatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matA.to_dgematrix() );
  
  for(CPPL_INT i=0; i<matB.m; i++){
    const CPPL_INT jmax =std::min(matB.n,i+matB.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax; j++){
      newmat(i,j)+=matB(i,j);
    }
  }
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}

//=============================================================================
/*! _dgsmatrix-_dgbmatrix operator */
inline _dgematrix operator-(const _dgsmatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matA.to_dgematrix() );
  
  for(CPPL_INT i=0; i<matB.m; i++){
    const CPPL_INT jmax =std::min(matB.n,i+matB.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax; j++){
      newmat(i,j)-=matB(i,j);
    }
  }
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}

//=============================================================================
/*! _dgsmatrix*_dgbmatrix operator */
inline _dgematrix operator*(const _dgsmatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matA.m, matB.n );
  newmat.zero();
  
  const std::vector<dcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<dcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    const CPPL_INT jmax =std::min(matB.n,it->j+matB.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),it->j-matB.kl); j<jmax; j++){
      newmat(it->i,j) += it->v*matB(it->j,j);
    }
  }
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}
