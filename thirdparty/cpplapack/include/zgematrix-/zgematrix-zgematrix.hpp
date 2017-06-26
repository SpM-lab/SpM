//=============================================================================
/*! zgematrix=zgematrix operator */
inline zgematrix& zgematrix::operator=(const zgematrix& mat)
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
/*! zgematrix+=zgematrix operator */
inline zgematrix& zgematrix::operator+=(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n || m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << m << "x" << n << ") += (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT size =m*n;
  for(CPPL_INT i=0; i<size; i++){
    array[i]+=mat.array[i];
  }
  return *this;
}

//=============================================================================
/*! zgematrix operator-= */
inline zgematrix& zgematrix::operator-=(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n || m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a sutraction." << std::endl
              << "Your input was (" << m << "x" << n << ") -= (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT size =m*n;
  for(CPPL_INT i=0; i<size; i++){
    array[i]-=mat.array[i];
  }
  return *this;
}

//=============================================================================
/*! zgematrix operator*= */
inline zgematrix& zgematrix::operator*=(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << m << "x" << n << ") *= (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat( m, mat.n );
  char transa ='n';
  char transb ='n';
  comple alpha =comple(1.,0.);
  comple beta =comple(0.,0.);
  
  zgemm_( &transa, &transb, &m, &mat.n, &n, &alpha, array, &m, mat.array, &mat.m, &beta, newmat.array, &m );
  
  swap(*this,newmat);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix+zgematrix operator */
inline _zgematrix operator+(const zgematrix& matA, const zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  zgematrix newmat(matA.m,matA.n);
  
  const CPPL_INT size =matA.m*matA.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =matA.array[i]+matB.array[i];
  }
  
  return _(newmat);
}

//=============================================================================
/*! zgematrix-zgematrix operator */
inline _zgematrix operator-(const zgematrix& matA, const zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  zgematrix newmat(matA.m,matA.n);
  
  const CPPL_INT size =matA.m*matA.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =matA.array[i]-matB.array[i];
  }
  
  return _(newmat);
}

//=============================================================================
/*! zgematrix*zgematrix operator */
inline _zgematrix operator*(const zgematrix& matA, const zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat( matA.m, matB.n );
  char transa ='n';
  char transb ='n';
  comple alpha =comple(1.,0.);
  comple beta =comple(0.,0.);
  
  zgemm_( &transa, &transb, &matA.m, &matB.n, &matA.n, &alpha, matA.array, &matA.m, matB.array, &matB.m, &beta, newmat.array, &matA.m );
  
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return Hadamerd product */
inline _zgematrix hadamerd(const zgematrix& matA, const zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( matA.m!=matB.m || matA.n!=matB.n ){
    ERROR_REPORT;
    std::cerr << "These two matrices can not make Hadamerd product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") and (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat(matA.m,matA.n);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =matA(i,j)*matB(i,j);
    }
  }
  
  return _(newmat);
}
