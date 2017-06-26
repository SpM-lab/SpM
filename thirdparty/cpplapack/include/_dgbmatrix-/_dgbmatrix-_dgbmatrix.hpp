//=============================================================================
/*! _dgbmatrix+_dgbmatrix operator */
inline _dgbmatrix operator+(const _dgbmatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if(matA.kl>matB.kl && matA.ku>matB.ku){
    for(CPPL_INT i=0; i<matB.m; i++){
      const CPPL_INT jmax =std::min(matB.n,i+matB.ku+1);
      for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax; j++){
        matA(i,j) += matB(i,j);
      }
    }
    
    matB.destroy();
    return matA;
  }

  else if(matB.kl>matA.kl && matB.ku>matA.ku){
    for(CPPL_INT i=0; i<matA.m; i++){
      const CPPL_INT jmax =std::min(matA.n,i+matA.ku+1);
      for(CPPL_INT j=std::max(CPPL_INT(0),i-matA.kl); j<jmax; j++){
        matB(i,j) += matA(i,j);
      }
    }
    
    matA.destroy();
    return matB;
  }
  
  else{
    dgbmatrix newmat(matA.m,matA.n,std::max(matA.kl,matB.kl),std::max(matA.ku,matB.ku));
    newmat.zero();
    
    for(CPPL_INT i=0; i<matA.m; i++){
      const CPPL_INT jmax1 =std::min(matA.n,i+matA.ku+1);
      for(CPPL_INT j=std::max(CPPL_INT(0),i-matA.kl); j<jmax1; j++){
        newmat(i,j) += matA(i,j);
      }
      const CPPL_INT jmax2 =std::min(matB.n,i+matB.ku+1);
      for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax2; j++){
        newmat(i,j) += matB(i,j);
      }
    }
    
    matA.destroy();
    matB.destroy();
    return _(newmat);
  }
}

//=============================================================================
/*! _dgbmatrix-_dgbmatrix operator */
inline _dgbmatrix operator-(const _dgbmatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if(matA.kl>matB.kl && matA.ku>matB.ku){
    for(CPPL_INT i=0; i<matB.m; i++){
      const CPPL_INT jmax =std::min(matB.n,i+matB.ku+1);
      for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax; j++){
        matA(i,j) -= matB(i,j);
      }
    }
    
    matB.destroy();
    return matA;
  }
  
  else{
    dgbmatrix newmat(matA.m,matA.n,std::max(matA.kl,matB.kl),std::max(matA.ku,matB.ku));
    newmat.zero();
    
    for(CPPL_INT i=0; i<matA.m; i++){
      const CPPL_INT jmax1 =std::min(matA.n,i+matA.ku+1);
      for(CPPL_INT j=std::max(CPPL_INT(0),i-matA.kl); j<jmax1; j++){
        newmat(i,j) += matA(i,j);
      }
      const CPPL_INT jmax2 =std::min(matB.n,i+matB.ku+1);
      for(CPPL_INT j=std::max(CPPL_INT(0),i-matB.kl); j<jmax2; j++){
        newmat(i,j) -= matB(i,j);
      }
    }
    
    matA.destroy();
    matB.destroy();
    return _(newmat);
  }
}

//=============================================================================
/*! _dgbmatrix*_dgbmatrix operator */
inline _dgbmatrix operator*(const _dgbmatrix& matA, const _dgbmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgbmatrix newmat( matA.m, matB.n, std::min(matA.kl+matB.kl,matA.m-1), std::min(matA.ku+matB.ku,matB.n-1) );
  newmat.zero();
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    const CPPL_INT jmax =std::min(newmat.n,i+newmat.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-newmat.kl); j<jmax; j++){
      const CPPL_INT kmax =std::min( std::min(matA.n,i+matA.ku+1), std::min(matB.m,j+matB.kl+1) );
      for(CPPL_INT k=std::max( std::max(CPPL_INT(0),i-matA.kl), std::max(CPPL_INT(0),j-matB.ku) ); k<kmax; k++){
        newmat(i,j) += matA(i,k)*matB(k,j);
      }
    }
  }
  
  matA.destroy();
  matB.destroy();
  return _(newmat);
}
