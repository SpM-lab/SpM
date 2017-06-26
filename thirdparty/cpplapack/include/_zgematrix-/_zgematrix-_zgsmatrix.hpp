//=============================================================================
/*! _zgematrix+_zgsmatrix operator */
inline _zgematrix operator+(const _zgematrix& matA, const _zgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.m!=matB.m || matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const size_t matB_data_size =matB.data.size();
  for(size_t c=0; c<matB_data_size; c++){
    const zcomponent& z =matB.data[c];
    matA(z.i,z.j) += z.v;
  }
  
  matB.destroy();
  return matA;
}

//=============================================================================
/*! _zgematrix-_zgsmatrix operator */
inline _zgematrix operator-(const _zgematrix& matA, const _zgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.m!=matB.m || matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const size_t matB_data_size =matB.data.size();
  for(size_t c=0; c<matB_data_size; c++){
    const zcomponent& z =matB.data[c];
    matA(z.i,z.j) -= z.v;
  }
  
  matB.destroy();
  return matA;
}

//=============================================================================
/*! _zgematrix*_zgsmatrix operator */
inline _zgematrix operator*(const _zgematrix& matA, const _zgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix newmat(matA.m, matB.n);
  newmat.zero();
  
  const size_t matB_data_size =matB.data.size();
  for(size_t c=0; c<matB_data_size; c++){
    const zcomponent& z =matB.data[c];
    for(CPPL_INT i=0; i<matA.m; i++){
      newmat(i,z.j) += matA(i,z.i)*z.v;
    }
  }
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}
