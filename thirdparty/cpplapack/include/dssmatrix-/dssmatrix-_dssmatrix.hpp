//=============================================================================
/*! dssmatrix=_dssmatrix operator */
inline dssmatrix& dssmatrix::operator=(const _dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
  shallow_copy(mat);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dssmatrix+=_dssmatrix operator */
inline dssmatrix& dssmatrix::operator+=(const _dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << n << "x" << n << ") += (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const std::vector<dcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    (*this)(it->i,it->j) +=it->v;
  }
  
  mat.destroy();
  return *this;
}

//=============================================================================
/*! dssmatrix-=_dssmatrix operator */
inline dssmatrix& dssmatrix::operator-=(const _dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a sutraction." << std::endl
              << "Your input was (" << n << "x" << n << ") -= (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const std::vector<dcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    (*this)(it->i,it->j) -=it->v;
  }
  
  mat.destroy();
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dssmatrix+_dssmatrix operator */
inline _dssmatrix operator+(const dssmatrix& matA, const _dssmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") + (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  dssmatrix newmat(matB);
  
  const std::vector<dcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<dcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    newmat(it->i,it->j) +=it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dssmatrix-_dssmatrix operator */
inline _dssmatrix operator-(const dssmatrix& matA, const _dssmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") - (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dssmatrix newmat(matB);
  
  const std::vector<dcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<dcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    newmat(it->i,it->j) -=it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dssmatrix*_dssmatrix operator */
/*
inline _dgsmatrix operator*(const dssmatrix& matA, const _dssmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.n << "x" << matA.n << ") * (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  dssmatrix newmat( matA.n, 0 );
  
  for(CPPL_INT c=0; c<matA.vol; c++){
    CPPL_INT k(matA.jndx[c]);
    std::vector< std::pair<CPPL_INT,CPPL_INT> >::iterator p;
    for(p=matB.line[k].begin(); p!=matB.line[k].end(); p++){
      newmat(matA.indx[c],p->first) +=matA.array[c]*matB.array[p->second];
    }
  }
  
  matB.destroy();
  return _(newmat);
}
*/
