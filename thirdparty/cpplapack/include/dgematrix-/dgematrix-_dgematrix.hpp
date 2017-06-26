//=============================================================================
/*! dgematrix=_dgematrix operator */
inline dgematrix& dgematrix::operator=(const _dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  shallow_copy(mat);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgematrix+=_dgematrix operator */
inline dgematrix& dgematrix::operator+=(const _dgematrix& mat)
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
  
  mat.destroy();
  return *this;
}

//=============================================================================
/*! dgematrix-=_dgematrix operator */
inline dgematrix& dgematrix::operator-=(const _dgematrix& mat)
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

  mat.destroy();
  return *this;
}

//=============================================================================
/*! dgematrix*=_dgematrix operator */
inline dgematrix& dgematrix::operator*=(const _dgematrix& mat)
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
  mat.destroy();
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgematrix+_dgematrix operator */
inline _dgematrix operator+(const dgematrix& matA, const _dgematrix& matB)
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
    matB.array[i] +=matA.array[i];
  }
  
  return matB;
}

//=============================================================================
/*! dgematrix-_dgematrix operator */
inline _dgematrix operator-(const dgematrix& matA, const _dgematrix& matB)
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
    matB.array[i] =matA.array[i]-matB.array[i];
  }
  
  return matB;
}

//=============================================================================
/*! dgematrix*_dgematrix operator */
inline _dgematrix operator*(const dgematrix& matA, const _dgematrix& matB)
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
  
  matB.destroy();
  return _(newmat);
}
