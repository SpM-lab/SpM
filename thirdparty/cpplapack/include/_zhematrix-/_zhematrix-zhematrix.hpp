//=============================================================================
/*! _zhematrix+zhematrix operator */
inline _zhematrix operator+(const _zhematrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT j=0; j<matA.n; j++){
    for(CPPL_INT i=j; i<matA.n; i++){
      matA.darray[j][i] += matB.darray[j][i];
    }
  }
  
  return matA;
}

//=============================================================================
/*! _zhematrix-zhematrix operator */
inline _zhematrix operator-(const _zhematrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") - (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT j=0; j<matA.n; j++){
    for(CPPL_INT i=j; i<matA.n; i++){
      matA.darray[j][i] -= matB.darray[j][i];
    }
  }
  
  return matA;
}

//=============================================================================
/*! _zhematrix*zhematrix operator */
inline _zgematrix operator*(const _zhematrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") * (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  matB.complete();
  
  zgematrix newmat( matA.n, matB.n );
  char side ='l';
  char uplo ='l';
  comple alpha =comple(1.,0.);
  comple beta =comple(0.,0.);
  
  zhemm_( &side, &uplo, &matA.n, &matB.n, &alpha, matA.array, &matA.n, matB.array, &matB.n, &beta, newmat.array, &newmat.m );
  
  matA.destroy();
  return _(newmat);
}
