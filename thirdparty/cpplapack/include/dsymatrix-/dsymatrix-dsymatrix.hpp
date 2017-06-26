//=============================================================================
/*! dsymatrix=dsymatrix operator */
inline dsymatrix& dsymatrix::operator=(const dsymatrix& mat)
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
/*! dsymatrix+=dsymatrix operator */
inline dsymatrix& dsymatrix::operator+=(const dsymatrix& mat)
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
      darray[j][i] +=mat.darray[j][i];
    }
  }
  
  return *this;
}

//=============================================================================
/*! dsymatrix operator-= */
inline dsymatrix& dsymatrix::operator-=(const dsymatrix& mat)
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
      darray[j][i] -=mat.darray[j][i];
    }
  }
  
  return *this;
}

//=============================================================================
/*! dsymatrix+dsymatrix operator */
inline _dsymatrix operator+(const dsymatrix& matA, const dsymatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dsymatrix newmat(matA.n);
  
  for(CPPL_INT j=0; j<matA.n; j++){
    for(CPPL_INT i=j; i<matA.n; i++){
      newmat.darray[j][i] =matA.darray[j][i] +matB.darray[j][i];
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! dsymatrix-dsymatrix operator */
inline _dsymatrix operator-(const dsymatrix& matA, const dsymatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") - (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dsymatrix newmat(matA.n);
  
  for(CPPL_INT j=0; j<matA.n; j++){
    for(CPPL_INT i=j; i<matA.n; i++){
      newmat.darray[j][i] =matA.darray[j][i] -matB.darray[j][i];
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! dsymatrix*dsymatrix operator */
inline _dgematrix operator*(const dsymatrix& matA, const dsymatrix& matB)
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
  
  dgematrix newmat(matA.n, matA.n);
  char side ='l';
  char uplo ='l';
  double alpha =1.;
  double beta =0.;
  
  dsymm_( &side, &uplo, &matA.n, &matB.n, &alpha, matA.array, &matA.n, matB.array, &matB.m, &beta, newmat.array, &newmat.m );
  
  return _(newmat);
}
