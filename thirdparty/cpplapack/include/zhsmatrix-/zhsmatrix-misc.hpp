//=============================================================================
/*! clear all the matrix data and set the sizes 0 */
inline void zhsmatrix::clear()
{CPPL_VERBOSE_REPORT;
  n =0;
  data.clear();
  line.clear();
}

//=============================================================================
/*! change the matrix into a zero matrix */
inline zhsmatrix& zhsmatrix::zero()
{CPPL_VERBOSE_REPORT;
  data.resize(0);
  for(CPPL_INT i=0; i<n; i++){ line[i].resize(0); }
  return *this;
}

//=============================================================================
/*! change the matrix into an identity matrix */
inline zhsmatrix& zhsmatrix::identity()
{CPPL_VERBOSE_REPORT;
  zero();
  for(CPPL_INT i=0; i<n; i++){
    put(i,i, comple(1.,0.));
  }
  return *this;
}

//=============================================================================
/*! change sign(+/-) of the matrix */
inline void zhsmatrix::chsign()
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator data_end =data.end();
  for(std::vector<zcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v =-it->v;
  }
}

//=============================================================================
/*! make a deep copy of the matrix */
inline void zhsmatrix::copy(const zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  n =mat.n;
  data =mat.data;
  line =mat.line;
}

//=============================================================================
/*! make a shallow copy of the matrix\n
  This function is not designed to be used in project codes. */
inline void zhsmatrix::shallow_copy(const _zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  data.clear();
  line.clear();
  
  n =mat.n;
  data.swap(mat.data);
  line.swap(mat.line);
  
  mat.nullify();
}

//=============================================================================
/*! resize the matrix */
inline zhsmatrix& zhsmatrix::resize(const CPPL_INT& _n, const CPPL_INT _c, const CPPL_INT _l)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _n<0 || _c<0 || _l<0 ){
    ERROR_REPORT;
    std::cerr << "Matrix sizes and the length of arrays must be positive integers. " << std::endl
              << "Your input was (" << _n << "," << _c << "," << _l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  n =_n;
  data.resize(0);
  data.reserve(_c);
  line.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    line[i].resize(0);
    line[i].reserve(_l);
  }
  
  return *this;
}

//=============================================================================
/*! stretch the matrix size */
inline void zhsmatrix::stretch(const CPPL_INT& dn)
{CPPL_VERBOSE_REPORT;
  if(dn==0){ return; }
  
#ifdef  CPPL_DEBUG
  if( n+dn<0 ){
    ERROR_REPORT;
    std::cerr << "The new matrix size must be larger than zero." << std::endl
              << "Your input was (" << dn << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  n +=dn;
  
  if(dn<0){
    //// delete components over the new size ////
    const std::vector<zcomponent>::reverse_iterator data_rend =data.rend();
    for(std::vector<zcomponent>::reverse_iterator it=data.rbegin(); it!=data_rend; it++){
      if( it->i>=n ){ del( CPPL_INT(data_rend-it-1) ); }
    }
    //// shrink line ////
    for(CPPL_INT i=0; i<-dn; i++){
      line.pop_back();
    }
  }
  else{//dn>0
    //// expand line ////
    for(CPPL_INT i=0; i<dn; i++){
      line.push_back( std::vector<CPPL_INT>(0) );
    }
  }
}

//=============================================================================
/*! check if the component is listed */
inline bool zhsmatrix::isListed(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT ii(std::max(i,j)), jj(std::min(i,j));
  
  const std::vector<CPPL_INT>::const_iterator line_ii_end =line[ii].end();
  for(std::vector<CPPL_INT>::const_iterator p=line[ii].begin(); p!=line_ii_end; p++){
    if(data[*p].j==jj){ return 1; }
  }
  
  return 0;
}

//=============================================================================
/*! return the element number of the component */
inline CPPL_INT zhsmatrix::number(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  const CPPL_INT ii(std::max(i,j)), jj(std::min(i,j));
  
  const std::vector<CPPL_INT>::const_iterator line_ii_end =line[ii].end();
  for(std::vector<CPPL_INT>::const_iterator p=line[ii].begin(); p!=line_ii_end; p++){
    if(data[*p].i==ii && data[*p].j==jj){
      return *p;
    }
  }
  
  return -1;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! get row of the matrix */
inline _zrovector zhsmatrix::row(const CPPL_INT& _m) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _m<0 || _m>m ){
    ERROR_REPORT;
    std::cerr << "Input row number must be between 0 and " << m << "." << std::endl
              << "Your input was " << _m << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zrovector vec( zrovector(n).zero() );
  
  const std::vector<CPPL_INT>::const_iterator line__m_end =line[_m].end();
  for(std::vector<CPPL_INT>::const_iterator p=line[_m].begin(); p!=line__m_end; p++){
    if( data[*p].i==_m ){//i>=j
      vec(data[*p].j) =data[*p].v;
    }
    else{//i<j
      vec(data[*p].i) =std::conj(data[*p].v);
    }
  }
  
  return _(vec);
}

//=============================================================================
/*! get column of the matrix */
inline _zcovector zhsmatrix::col(const CPPL_INT& _n) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _n<0 || _n>n ){
    ERROR_REPORT;
    std::cerr << "Input row number must be between 0 and " << n << "." << std::endl
              << "Your input was " << _n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zcovector vec( zcovector(m).zero() );
  
  const std::vector<CPPL_INT>::const_iterator line__n_end =line[_n].end();
  for(std::vector<CPPL_INT>::const_iterator p=line[_n].begin(); p!=line__n_end; p++){
    if( data[*p].i==_n ){//i<j
      vec(data[*p].j) =std::conj(data[*p].v);
    }
    else{//i>=j
      vec(data[*p].i) =data[*p].v;
    }
  }
  
  return _(vec);
}

//=============================================================================
/*! erase components less than DBL_MIN */
inline void zhsmatrix::diet(const double eps)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::reverse_iterator data_rend =data.rend();
  for(std::vector<zcomponent>::reverse_iterator it=data.rbegin(); it!=data_rend; it++){
    if( fabs(it->v.real())<eps && fabs(it->v.imag())<eps ){
      del( CPPL_INT(data_rend-it-1) );
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! swap two matrices */
inline void swap(zhsmatrix& A, zhsmatrix& B)
{CPPL_VERBOSE_REPORT;
  std::swap(A.n,B.n);
  std::swap(A.data,B.data);
  std::swap(A.line,B.line);
}

//=============================================================================
/*! convert user object to smart-temporary object */
inline _zhsmatrix _(zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  _zhsmatrix newmat;
  //////// shallow copy ////////
  newmat.n =mat.n;
  std::swap(newmat.data, mat.data);
  std::swap(newmat.line, mat.line);
  //////// nullify ////////
  mat.n =0;
  return newmat;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! health checkup */
inline void zhsmatrix::checkup()
{CPPL_VERBOSE_REPORT;
  //////// complex diagonal ////////
  for(CPPL_INT i=0; i<m; i++){
    if( std::fabs((*this)(i,i).imag()) > DBL_MIN ){
      ERROR_REPORT;
      std::cerr << "Diagonal components of a Hermitian matrix have to be real numbers." << std::endl
                << "(*this)(" << i << "," << i << ") was a complex number, " << (*this)(i,i) << "." << std::endl;
      exit(1);
    }
  }
  
  //////// NOTE ////////
  std::cerr << "# [NOTE]@zhsmatrix::checkup(): This symmetric sparse matrix is fine." << std::endl;
}
