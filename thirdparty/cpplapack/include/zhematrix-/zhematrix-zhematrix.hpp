//=============================================================================
/*! zhematrix=zhematrix operator */
inline zhematrix& zhematrix::operator=(const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  if(array!=mat.array){ // if it is NOT self substitution
    copy(mat);
  }
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix+=zhematrix operator */
inline zhematrix& zhematrix::operator+=(const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << n << "x" << n << ") += (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      darray[j][i] += mat.darray[j][i];
    }
  }
  
  return *this;
}

//=============================================================================
/*! zhematrix operator-= */
inline zhematrix& zhematrix::operator-=(const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a sutraction." << std::endl
              << "Your input was (" << n << "x" << n << ") -= (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      darray[j][i] -= mat.darray[j][i];
    }
  }
  
  return *this;
}

//=============================================================================
/*! zhematrix+zhematrix operator */
inline _zhematrix operator+(const zhematrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zhematrix newmat(matA.n);
  
  for(CPPL_INT j=0; j<matA.n; j++){
    for(CPPL_INT i=j; i<matA.n; i++){
      newmat.darray[j][i] =matA.darray[j][i] +matB.darray[j][i];
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! zhematrix-zhematrix operator */
inline _zhematrix operator-(const zhematrix& matA, const zhematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") - (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zhematrix newmat(matA.n);
  
  for(CPPL_INT j=0; j<matA.n; j++){
    for(CPPL_INT i=j; i<matA.n; i++){
      newmat.darray[j][i] =matA.darray[j][i] -matB.darray[j][i];
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! zhematrix*zhematrix operator */
inline _zgematrix operator*(const zhematrix& matA, const zhematrix& matB)
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
  
  return _(newmat);
}
