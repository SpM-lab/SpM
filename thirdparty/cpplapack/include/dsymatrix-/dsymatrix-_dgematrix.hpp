//=============================================================================
/*! _dgematrix+dsymatrix operator */
inline _dgematrix operator+(const dsymatrix& matA, const _dgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<matA.n; i++) {
    for(CPPL_INT j=0; j<matA.n; j++) {
      matB(i,j) += matA(i,j);
    }
  }
  
  return matB;
}

//=============================================================================
/*! _dgematrix-dgematrix operator */
inline _dgematrix operator-(const dsymatrix& matA, const _dgematrix& matB)
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
  
  return matB;
}

//=============================================================================
/*! _dgematrix*dgematrix operator */
inline _dgematrix operator*(const dsymatrix& matA, const _dgematrix& matB)
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
  char side ='l';
  char uplo ='l';
  double alpha =1.;
  double beta =0.;

  dsymm_( &side, &uplo, &matA.n, &matB.n, &alpha, matA.array, &matA.n, matB.array, &matB.m, &beta, newmat.array, &newmat.m );
  
  matB.destroy();
  return _(newmat);
}
