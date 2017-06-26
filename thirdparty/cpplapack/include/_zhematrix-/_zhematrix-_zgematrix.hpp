//=============================================================================
/*! _zgematrix+zhematrix operator */
inline _zgematrix operator+(const _zhematrix& matA, const _zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<matA.n; i++){
    for(CPPL_INT j=0; j<matA.n; j++){
      matB(i,j) += matA(i,j);
    }
  }
  
  matA.destroy();
  return matB;
}

//=============================================================================
/*! _zhematrix-zgematrix operator */
inline _zgematrix operator-(const _zhematrix& matA, const _zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<matA.n; i++){
    for(CPPL_INT j=0; j<matA.n; j++){
      matB(i,j) =matA(i,j)-matB(i,j);
    }
  }
  
  matA.destroy();
  return matB;
}

//=============================================================================
/*! _zgematrix*zgematrix operator */
inline _zgematrix operator*(const _zhematrix& matA, const _zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat( matA.n, matB.n );
  char side ='l';
  char uplo ='l';
  comple alpha =comple(1.,0.);
  comple beta =comple(0.,0.);
  
  zhemm_( &side, &uplo, &matA.n, &matB.n, &alpha, matA.array, &matA.n, matB.array, &matB.m, &beta, newmat.array, &newmat.m );
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}
