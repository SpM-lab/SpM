//=============================================================================
/*! dsymatrix=_dsymatrix operator */
inline dsymatrix& dsymatrix::operator=(const _dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  shallow_copy(mat);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dsymatrix+=_dsymatrix operator */
inline dsymatrix& dsymatrix::operator+=(const _dsymatrix& mat)
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
  
  mat.destroy();
  return *this;
}

//=============================================================================
/*! dsymatrix-=_dsymatrix operator */
inline dsymatrix& dsymatrix::operator-=(const _dsymatrix& mat)
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
  
  mat.destroy();
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dsymatrix+_dsymatrix operator */
inline _dsymatrix operator+(const dsymatrix& matA, const _dsymatrix& matB)
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
      matB.darray[j][i] +=matA.darray[j][i];
    }
  }
  
  return matB;
}

//=============================================================================
/*! dsymatrix-_dsymatrix operator */
inline _dsymatrix operator-(const dsymatrix& matA, const _dsymatrix& matB)
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
      matB.darray[j][i] =matA.darray[j][i] -matB.darray[j][i];
    }
  }
  
  return matB;
}

//=============================================================================
/*! dsymatrix*_dsymatrix operator */
inline _dgematrix operator*(const dsymatrix& matA, const _dsymatrix& matB)
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
  
  dgematrix newmat( matA.n, matA.n );
  char side ='l';
  char uplo ='l';
  double alpha =1.;
  double beta =0.;
  
  dsymm_( &side, &uplo, &matA.n, &matB.n, &alpha, matA.array, &matA.n, matB.array, &matB.m, &beta, newmat.array, &newmat.m );
  
  matB.destroy();
  return _(newmat);
}
