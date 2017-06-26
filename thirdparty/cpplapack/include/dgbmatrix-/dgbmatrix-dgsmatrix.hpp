//=============================================================================
/*! dgbmatrix+dgsmatrix operator */
inline _dgematrix operator+(const dgbmatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matB.to_dgematrix() );
  
  for(CPPL_INT i=0; i<matA.m; i++){
    const CPPL_INT jmax =std::min(matA.n,i+matA.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matA.kl); j<jmax; j++){
      newmat(i,j)+=matA(i,j);
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgbmatrix-dgsmatrix operator */
inline _dgematrix operator-(const dgbmatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( -matB.to_dgematrix() );
  
  for(CPPL_INT i=0; i<matA.m; i++){
    const CPPL_INT jmax =std::min(matA.n,i+matA.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matA.kl); j<jmax; j++){
      newmat(i,j)-=matA(i,j);
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgbmatrix*dgsmatrix operator */
inline _dgematrix operator*(const dgbmatrix& matA, const dgsmatrix& matB)
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
  
  const std::vector<dcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<dcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    const CPPL_INT imax =std::min(matA.m,it->i+matA.kl);
    for(CPPL_INT i=std::max(CPPL_INT(0),it->i-(matA.ku+1)); i<imax; i++){
      newmat(i,it->j) += matA(i,it->i)*it->v;
    }
  }
  
  return _(newmat);
}
