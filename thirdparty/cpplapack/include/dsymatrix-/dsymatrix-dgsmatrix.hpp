//=============================================================================
/*! dsymatrix+dgsmatrix operator */
inline _dgematrix operator+(const dsymatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m || matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matA.to_dgematrix() );
  
  const std::vector<dcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<dcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    newmat(it->i,it->j) += it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dsymatrix-dgsmatrix operator */
inline _dgematrix operator-(const dsymatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m || matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat( matA.to_dgematrix() );
  
  const std::vector<dcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<dcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    newmat(it->i,it->j) -= it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dsymatrix*dgsmatrix operator */
inline _dgematrix operator*(const dsymatrix& matA, const dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix newmat(matA.n, matB.n);
  newmat.zero();
  
  const std::vector<dcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<dcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    for(CPPL_INT i=0; i<matA.n; i++){
      newmat(i,it->j) += matB(i,it->i)*it->v;
    }
  }
  
  return _(newmat);
}
