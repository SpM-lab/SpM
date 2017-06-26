//=============================================================================
/*! dgematrix=dgematrix operator */
inline dgematrix& dgematrix::operator=(const dgematrix& mat)
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
/*! dgematrix+=dgematrix operator */
inline dgematrix& dgematrix::operator+=(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n || m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << m << "x" << n << ") += (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT mn =m*n;
  for(CPPL_INT i=0; i<mn; i++){
    array[i]+=mat.array[i];
  }
  return *this;
}

//=============================================================================
/*! dgematrix operator-= */
inline dgematrix& dgematrix::operator-=(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n || m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a sutraction." << std::endl
              << "Your input was (" << m << "x" << n << ") -= (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT mn =m*n;
  for(CPPL_INT i=0; i<mn; i++){
    array[i]-=mat.array[i];
  }
  return *this;
}

//=============================================================================
/*! dgematrix operator*= */
inline dgematrix& dgematrix::operator*=(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << m << "x" << n << ") *= (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( m, mat.n );
  char transa ='n';
  char transb ='n';
  double alpha =1.;
  double beta =0.;
  
  dgemm_( &transa, &transb, &m, &mat.n, &n, &alpha, array, &m, mat.array, &mat.m, &beta, newmat.array, &m );
  
  swap(*this,newmat);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgematrix+dgematrix operator */
inline _dgematrix operator+(const dgematrix& matA, const dgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  dgematrix newmat(matA.m,matA.n);
  
  const CPPL_INT mn =newmat.m*newmat.n;
  for(CPPL_INT i=0; i<mn; i++){
    newmat.array[i] =matA.array[i]+matB.array[i];
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgematrix-dgematrix operator */
inline _dgematrix operator-(const dgematrix& matA, const dgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  dgematrix newmat(matA.m,matA.n);
  
  const CPPL_INT mn =newmat.m*newmat.n;
  for(CPPL_INT i=0; i<mn; i++){
    newmat.array[i] =matA.array[i]-matB.array[i];
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgematrix*dgematrix operator */
inline _dgematrix operator*(const dgematrix& matA, const dgematrix& matB)
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
  
  return _(newmat);
}

//=============================================================================
/*! dgematrix%dgematrix operator */
inline _drovector operator%(const dgematrix& matA, const dgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.m!=matB.m || matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") % (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  drovector newvec( matA.n );
  
  newvec.zero();
  for(CPPL_INT j=0; j<matA.n; j++){
    for(CPPL_INT i=0; i<matA.m; i++){
      newvec(j) +=matA(i,j)*matB(i,j);
    }
  }
  
  return _(newvec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return Hadamerd product */
inline _dgematrix hadamerd(const dgematrix& matA, const dgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( matA.m!=matB.m || matA.n!=matB.n ){
    ERROR_REPORT;
    std::cerr << "These two matrices can not make Hadamerd product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") and (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  dgematrix newmat(matA.m,matA.n);
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =matA(i,j)*matB(i,j);
    }
  }
  return _(newmat);
}
