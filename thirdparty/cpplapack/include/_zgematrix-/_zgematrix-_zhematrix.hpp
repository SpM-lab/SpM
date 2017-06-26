//=============================================================================
/*! _zgematrix+_zhematrix operator */
inline _zgematrix operator+(const _zgematrix& matA, const _zhematrix& matB)
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
      matA(i,j) += matB(i,j);
    }
  }
  
  matB.destroy();
  return matA;
}

//=============================================================================
/*! _zgematrix-_zhematrix operator */
inline _zgematrix operator-(const _zgematrix& matA, const _zhematrix& matB)
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
      matA(i,j) -= matB(i,j);
    }
  }
  
  matB.destroy();
  return matA;
}

//=============================================================================
/*! _zgematrix*_zhematrix operator */
inline _zgematrix operator*(const _zgematrix& matA, const _zhematrix& matB)
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
  char side ='R';
  char uplo ='l';
  comple alpha =comple(1.,0.);
  comple beta =comple(0.,0.);

  zhemm_( &side, &uplo, &newmat.m, &newmat.n, &alpha, matB.array, &matB.n, matA.array, &matA.m, &beta, newmat.array, &newmat.m );
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}
