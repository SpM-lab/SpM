//=============================================================================
/*! _dsymatrix+_dgbmatrix operator */
inline _dgematrix operator+(const _dsymatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat(matA.n, matA.n);
  
  for(CPPL_INT i=0; i<matB.m; i++){
    for(CPPL_INT j=i; j<matA.n; j++){
      newmat(i,j) = newmat(j,i) = matA(i,j);
    }
    const CPPL_INT jmax =std::min(matB.n,i+matB.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax; j++){
      newmat(i,j) += matB(i,j);
    }
  }
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}

//=============================================================================
/*! _dsymatrix-_dgbmatrix operator */
inline _dgematrix operator-(const _dsymatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat(matA.n, matA.n);
  
  for(CPPL_INT i=0; i<matB.m; i++){
    for(CPPL_INT j=i; j<matA.n; j++){
      newmat(i,j) = newmat(j,i) = matA(i,j);
    }
    const CPPL_INT jmax =std::min(matB.n,i+matB.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax; j++){
      newmat(i,j) -= matB(i,j);
    }
  }

  matA.destroy();
  matB.destroy();
  return _(newmat);
}

//=============================================================================
/*! _dsymatrix*_dgbmatrix operator */
inline _dgematrix operator*(const _dsymatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matA.n, matB.n );
  newmat.zero();
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      const CPPL_INT kmax =std::min(matB.m,j+matB.kl+1);
      for(CPPL_INT k=std::max(CPPL_INT(0),j-matB.ku); k<kmax; k++){
        newmat(i,j) += matA(i,k)*matB(k,j);
      }
    }
  }
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}
