//=============================================================================
/*! zgematrix=_zgematrix operator */
inline zgematrix& zgematrix::operator=(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  shallow_copy(mat);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix+=_zgematrix operator */
inline zgematrix& zgematrix::operator+=(const _zgematrix& mat)
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
    array[i] += mat.array[i];
  }
  
  mat.destroy();
  return *this;
}

//=============================================================================
/*! zgematrix-=_zgematrix operator */
inline zgematrix& zgematrix::operator-=(const _zgematrix& mat)
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
    array[i] -= mat.array[i];
  }

  mat.destroy();
  return *this;
}

//=============================================================================
/*! zgematrix*=_zgematrix operator */
inline zgematrix& zgematrix::operator*=(const _zgematrix& mat)
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
  mat.destroy();
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix+_zgematrix operator */
inline _zgematrix operator+(const zgematrix& matA, const _zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT size =matA.m*matA.n;
  for(CPPL_INT i=0; i<size; i++){
    matB.array[i] +=matA.array[i];
  }
  
  return matB;
}

//=============================================================================
/*! zgematrix-_zgematrix operator */
inline _zgematrix operator-(const zgematrix& matA, const _zgematrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  const CPPL_INT size =matA.m*matA.n;
  for(CPPL_INT i=0; i<size; i++){
    matB.array[i] =matA.array[i]-matB.array[i];
  }
  
  return matB;
}

//=============================================================================
/*! zgematrix*_zgematrix operator */
inline _zgematrix operator*(const zgematrix& matA, const _zgematrix& matB)
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
  
  matB.destroy();
  return _(newmat);
}
