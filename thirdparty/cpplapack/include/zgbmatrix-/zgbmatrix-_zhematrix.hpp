//=============================================================================
/*! zgbmatrix+_zhematrix operator */
inline _zgematrix operator+(const zgbmatrix& matA, const _zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat(matB.n,matB.n);
  
  for(CPPL_INT i=0; i<matA.m; i++){
    for(CPPL_INT j=0; j<matB.n; j++){
      newmat(i,j) =matB(i,j);
    }
    const CPPL_INT jmax =std::min(matA.n,i+matA.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matA.kl); j<jmax; j++){
      newmat(i,j) += matA(i,j);
    }
  }
  
  matB.destroy();
  return _(newmat);
}

//=============================================================================
/*! zgbmatrix-_zhematrix operator */
inline _zgematrix operator-(const zgbmatrix& matA, const _zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat(matB.n,matB.n);
  
  for(CPPL_INT i=0; i<matA.m; i++){
    for(CPPL_INT j=0; j<matB.n; j++){
      newmat(i,j) =-matB(i,j);
    }
    const CPPL_INT jmax =std::min(matA.n,i+matA.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matA.kl); j<jmax; j++){
      newmat(i,j) +=matA(i,j);
    }
  }
  
  matB.destroy();
  return _(newmat);
}

//=============================================================================
/*! zgbmatrix*_zhematrix operator */
inline _zgematrix operator*(const zgbmatrix& matA, const _zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat( matA.m, matB.n );
  newmat.zero();
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      const CPPL_INT kmax =std::min(matA.n,i+matA.ku+1);
      for(CPPL_INT k=std::max(CPPL_INT(0),i-matA.kl); k<kmax; k++){
        newmat(i,j)+=matA(i,k)*matB(k,j);
      }
    }
  }
  
  matB.destroy();
  return _(newmat);
}
