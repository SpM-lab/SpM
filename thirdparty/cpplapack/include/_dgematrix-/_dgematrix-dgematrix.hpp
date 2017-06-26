//=============================================================================
/*! _dgematrix+dgematrix operator */
inline _dgematrix operator+(const _dgematrix& matA, const dgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT mn =matA.m*matA.n;
  for(CPPL_INT i=0; i<mn; i++){
    matA.array[i]+=matB.array[i];
  }
  
  return matA;
}

//=============================================================================
/*! _dgematrix-dgematrix operator */
inline _dgematrix operator-(const _dgematrix& matA, const dgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  const CPPL_INT mn =matA.m*matA.n;
  for(CPPL_INT i=0; i<mn; i++){
    matA.array[i]-=matB.array[i];
  }
  
  return matA;
}

//=============================================================================
/*! _dgematrix*dgematrix operator */
inline _dgematrix operator*(const _dgematrix& matA, const dgematrix& matB)
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
  char transa ='n';
  char transb ='n';
  double alpha =1.;
  double beta =0.;
  
  dgemm_( &transa, &transb, &matA.m, &matB.n, &matA.n, &alpha, matA.array, &matA.m, matB.array, &matB.m, &beta, newmat.array, &matA.m );
  
  matA.destroy();
  return _(newmat);
}
