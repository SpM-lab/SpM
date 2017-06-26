//=============================================================================
/*! _dgematrix+dsymatrix operator */
inline _dgematrix operator+(const _dgematrix& matA, const dsymatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<matB.n; i++){
    for(CPPL_INT j=0; j<matB.n; j++){
      matA(i,j)+=matB(i,j);
    }
  }
  
  return matA;
}

//=============================================================================
/*! _dgematrix-dsymatrix operator */
inline _dgematrix operator-(const _dgematrix& matA, const dsymatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<matB.n; i++){
    for(CPPL_INT j=0; j<matB.n; j++){
      matA(i,j)-=matB(i,j);
    }
  }
  
  return matA;
}

//=============================================================================
/*! _dgematrix*dsymatrix operator */
inline _dgematrix operator*(const _dgematrix& matA, const dsymatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matA.m, matB.n );
  char side ='R';
  char uplo ='l';
  double alpha =1.;
  double beta =0.;
  
  dsymm_( &side, &uplo, &newmat.m, &newmat.n, &alpha, matB.array, &matB.n, matA.array, &matA.m, &beta, newmat.array, &newmat.m );
  
  matA.destroy();
  return _(newmat);
}
