//=============================================================================
/*! zhsmatrix=zhsmatrix operator */
inline zhsmatrix& zhsmatrix::operator=(const zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  if(&data!=&mat.data){ // if it is NOT self substitution
    copy(mat);
  }
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhsmatrix+=zhsmatrix operator */
inline zhsmatrix& zhsmatrix::operator+=(const zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << n << "x" << n << ") += (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const std::vector<zcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    (*this)(it->i,it->j) +=it->v;
  }
  
  return *this;
}

//=============================================================================
/*! zhsmatrix-=zhsmatrix operator */
inline zhsmatrix& zhsmatrix::operator-=(const zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a sutraction." << std::endl
              << "Your input was (" << n << "x" << n << ") -= (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const std::vector<zcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<zcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    (*this)(it->i,it->j) -=it->v;
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhsmatrix+zhsmatrix operator */
inline _zhsmatrix operator+(const zhsmatrix& matA, const zhsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zhsmatrix newmat(matA);
  
  const std::vector<zcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<zcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    newmat(it->i,it->j) +=it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zhsmatrix-zhsmatrix operator */
inline _zhsmatrix operator-(const zhsmatrix& matA, const zhsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") - (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zhsmatrix newmat(matA);
  
  const std::vector<zcomponent>::const_iterator matB_data_end =matB.data.end();
  for(std::vector<zcomponent>::const_iterator it=matB.data.begin(); it!=matB_data_end; it++){
    newmat(it->i,it->j) -= it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zhsmatrix*zhsmatrix operator */
/*
inline _zgsmatrix operator*(const zhsmatrix& matA, const zhsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") * (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zhsmatrix newmat( matA.n, 0 );
  
  for(CPPL_INT c=0; c<matA.vol; c++){
    CPPL_INT k(matA.jndx[c]);
    std::vector< std::pair<CPPL_INT,CPPL_INT> >::iterator p;
    for(p=matB.Col[k].begin(); p!=matB.Col[k].end(); p++){
      newmat(matA.indx[c],p->first) +=matA.array[c]*matB.array[p->second];
    }
  }
  
  return _(newmat);
}
*/
