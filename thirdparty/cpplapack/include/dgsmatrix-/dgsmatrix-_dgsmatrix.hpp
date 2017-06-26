//=============================================================================
/*! dgsmatrix=_dgsmatrix operator */
inline dgsmatrix& dgsmatrix::operator=(const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  shallow_copy(mat);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgsmatrix+=_dgsmatrix operator */
inline dgsmatrix& dgsmatrix::operator+=(const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n || m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << m << "x" << n << ") += (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const std::vector<dcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    (*this)(it->i,it->j) += it->v;
  }
  
  mat.destroy();
  return *this;
}

//=============================================================================
/*! dgsmatrix-=_dgsmatrix operator */
inline dgsmatrix& dgsmatrix::operator-=(const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n || m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a sutraction." << std::endl
              << "Your input was (" << m << "x" << n << ") -= (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const std::vector<dcomponent>::const_iterator mat_data_end =mat.data.end();
  for(std::vector<dcomponent>::const_iterator it=mat.data.begin(); it!=mat_data_end; it++){
    (*this)(it->i,it->j) -= it->v;
  }
  
  mat.destroy();
  return *this;
}

//=============================================================================
/*! dgsmatrix*=_dgsmatrix operator */
inline dgsmatrix& dgsmatrix::operator*=(const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << m << "x" << n << ") *= (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgsmatrix newmat(m, mat.n);
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    CPPL_INT k =it->j;
    const std::vector<CPPL_INT>::const_iterator mat_rows_k_end =mat.rows[k].end();
    for(std::vector<CPPL_INT>::const_iterator p=mat.rows[k].begin(); p!=mat_rows_k_end; p++){
      newmat(it->i,mat.data[*p].j) += it->v*mat.data[*p].v;
    }
  }
  
  swap(*this,newmat);
  mat.destroy();
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgsmatrix+_dgsmatrix operator */
inline _dgsmatrix operator+(const dgsmatrix& matA, const _dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a summation." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") + (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  dgsmatrix newmat(matB);
  
  const std::vector<dcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<dcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    newmat(it->i,it->j) += it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgsmatrix-_dgsmatrix operator */
inline _dgsmatrix operator-(const dgsmatrix& matA, const _dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.n || matA.m!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a subtraction." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") - (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgsmatrix newmat(matB);
  newmat.chsign();
  
  const std::vector<dcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<dcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    newmat(it->i,it->j) += it->v;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgsmatrix*_dgsmatrix operator */
inline _dgsmatrix operator*(const dgsmatrix& matA, const _dgsmatrix& matB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matA.n!=matB.m){
    ERROR_REPORT;
    std::cerr << "These two matrises can not make a product." << std::endl
              << "Your input was (" << matA.m << "x" << matA.n << ") * (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgsmatrix newmat(matA.m, matB.n);

  const std::vector<dcomponent>::const_iterator matA_data_end =matA.data.end();
  for(std::vector<dcomponent>::const_iterator it=matA.data.begin(); it!=matA_data_end; it++){
    CPPL_INT k =it->j;
    const std::vector<CPPL_INT>::iterator matB_rows_k_end =matB.rows[k].end();
    for(std::vector<CPPL_INT>::iterator p=matB.rows[k].begin(); p!=matB_rows_k_end; p++){
      newmat(it->i,matB.data[*p].j) += it->v*matB.data[*p].v;
    }
  }
  
  matB.destroy();
  return _(newmat);
}
