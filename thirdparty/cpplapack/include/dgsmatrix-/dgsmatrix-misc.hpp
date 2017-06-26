//=============================================================================
/*! clear all the matrix data and set the sizes 0 */
inline void dgsmatrix::clear()
{CPPL_VERBOSE_REPORT;
  m =0;
  n =0;
  data.clear();
  rows.clear();
  cols.clear();
}

//=============================================================================
/*! change the matrix into a zero matrix */
inline dgsmatrix& dgsmatrix::zero()
{CPPL_VERBOSE_REPORT;
  data.resize(0);
  for(CPPL_INT i=0; i<m; i++){ rows[i].resize(0); }
  for(CPPL_INT j=0; j<n; j++){ cols[j].resize(0); }
  return *this;
}

//=============================================================================
/*! change the matrix into an identity matrix */
inline dgsmatrix& dgsmatrix::identity()
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n){
    ERROR_REPORT;
    std::cerr << "Only square matrix can be a identity matrix." << std::endl
              << "The matrix size was " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zero();
  for(CPPL_INT i=0; i<m; i++){
    put(i,i,1.);
  }
  return *this;
}

//=============================================================================
/*! change sign(+/-) of the matrix */
inline void dgsmatrix::chsign()
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator data_end =data.end();
  for(std::vector<dcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v =-it->v;
  }
}

//=============================================================================
/*! make a deep copy of the matrix */
inline void dgsmatrix::copy(const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  m =mat.m;
  n =mat.n;
  data =mat.data;
  rows =mat.rows;
  cols =mat.cols;
}

//=============================================================================
/*! make a shallow copy of the matrix\n
  This function is not designed to be used in project codes. */
inline void dgsmatrix::shallow_copy(const _dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  data.clear();
  rows.clear();
  cols.clear();
  
  m =mat.m;
  n =mat.n;
  data.swap(mat.data);
  rows.swap(mat.rows);
  cols.swap(mat.cols);
  
  mat.nullify();
}

//=============================================================================
/*! resize the matrix */
inline dgsmatrix& dgsmatrix::resize(const CPPL_INT& _m, const CPPL_INT& _n, const CPPL_INT _c, const CPPL_INT _l)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _m<0 || _n<0 || _c<0 ){
    ERROR_REPORT;
    std::cerr << "Matrix sizes and the length of arrays must be positive integers. " << std::endl
              << "Your input was (" << _m << "," << _n << "," << _c << "," << _l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  m =_m;
  n =_n;
  data.resize(0);
  data.reserve(_c);
  rows.resize(m);
  for(CPPL_INT i=0; i<m; i++){
    rows[i].resize(0);
    rows[i].reserve(_l);
  }
  cols.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    cols[i].resize(0);
    cols[i].reserve(_l);
  }
  
  return *this;
}

//=============================================================================
/*! stretch the matrix size */
inline void dgsmatrix::stretch(const CPPL_INT& dm, const CPPL_INT& dn)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( m+dm<0 || n+dn<0 ){
    ERROR_REPORT;
    std::cerr << "The new matrix size must be larger than zero. " << std::endl
              << "Your input was (" << dm << ", " << dn << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// zero ////////
  if(dm==0 && dn==0){ return; }
  
  //////// non-zero ////////
  m +=dm;
  n +=dn;
  
  //// for rows ////
  if(dm<0){
    //// delete components over the new size ////
    const std::vector<dcomponent>::reverse_iterator data_rend =data.rend();
    for(std::vector<dcomponent>::reverse_iterator it=data.rbegin(); it!=data_rend; it++){
      if( it->i>=m ){ del( CPPL_INT(data_rend-it-1) ); }
    }
    //// shrink rows ////
    for(CPPL_INT i=0; i<-dm; i++){
      rows.pop_back();
    }
  }
  else{//dm>0
    //// expand rows ////
    for(CPPL_INT i=0; i<dm; i++){
      rows.push_back( std::vector<CPPL_INT>(0) );
    }
  }
  
  //// for cols ////
  if(dn<0){
    //// delete components over the new size ////
    const std::vector<dcomponent>::reverse_iterator data_rend =data.rend();
    for(std::vector<dcomponent>::reverse_iterator it=data.rbegin(); it!=data_rend; it++){
      if( it->j>=n ){ del( CPPL_INT(data_rend-it-1) ); }
    }
    for(CPPL_INT j=0; j<-dn; j++){
      cols.pop_back();
    }
  }
  else{//dn>0
    //// expand cols ////
    for(CPPL_INT j=0; j<dn; j++){
      cols.push_back( std::vector<CPPL_INT>(0) );
    }
  }
}

//=============================================================================
/*! check if the component is listed */
inline bool dgsmatrix::isListed(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  const std::vector<CPPL_INT>::const_iterator rows_i_end =rows[i].end();
  for(std::vector<CPPL_INT>::const_iterator p=rows[i].begin(); p!=rows_i_end; p++){
    if(data[*p].j==j){ return 1; }
  }
  
  return 0;
}


//=============================================================================
/*! return the element number of the component */
inline CPPL_INT dgsmatrix::number(const CPPL_INT& i, const CPPL_INT& j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  const std::vector<CPPL_INT>::iterator rows_i_end =rows[i].end();
  for(std::vector<CPPL_INT>::iterator p=rows[i].begin(); p!=rows_i_end; p++){
    if(data[*p].j==j){ return *p; }
  }
  
  return -1;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! get row of the matrix */
inline _drovector dgsmatrix::row(const CPPL_INT& _m) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _m<0 || _m>m ){
    ERROR_REPORT;
    std::cerr << "Input row number must be between 0 and " << m << "." << std::endl
              << "Your input was " << _m << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  drovector vec(n);
  vec.zero();
  
  const std::vector<CPPL_INT>::const_iterator rows__m_end =rows[_m].end();
  for(std::vector<CPPL_INT>::const_iterator p=rows[_m].begin(); p!=rows__m_end; p++){
    vec(data[*p].j) =data[*p].v;
  }
  
  return _(vec);
}

//=============================================================================
/*! get column of the matrix */
inline _dcovector dgsmatrix::col(const CPPL_INT& _n) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _n<0 || _n>n ){
    ERROR_REPORT;
    std::cerr << "Input row number must be between 0 and " << n << "." << std::endl
              << "Your input was " << _n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dcovector vec(m);
  vec.zero();
  
  const std::vector<CPPL_INT>::const_iterator cols__n_end =cols[_n].end();
  for(std::vector<CPPL_INT>::const_iterator p=cols[_n].begin(); p!=cols__n_end; p++){
    vec(data[*p].i) =data[*p].v;
  }
  
  return _(vec);
}

//=============================================================================
/*! erase components less than DBL_MIN */
inline void dgsmatrix::diet(const double eps)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::reverse_iterator data_rend =data.rend();
  for(std::vector<dcomponent>::reverse_iterator it=data.rbegin(); it!=data_rend; it++){
    if( fabs(it->v)<eps ){ del( CPPL_INT(data_rend-it-1) ); }
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! health checkup */
inline void dgsmatrix::checkup()
{CPPL_VERBOSE_REPORT;
  //////// write ////////
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    std::cerr << "array[" << it-data.begin() << "] = (" << it->i << "," << it->j << ") = " << it->v << std::endl;
  }
  std::cerr << std::endl;
  
  for(CPPL_INT i=0; i<m; i++){
    std::cerr << "rows[" << i << "] =" << std::flush;
    const size_t rows_i_size =rows[i].size();
    for(size_t k=0; k<rows_i_size; k++){
      std::cerr << " " << rows[i][k] << std::flush;
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  
  for(CPPL_INT j=0; j<n; j++){
    std::cerr << "cols[" << j << "] =" << std::flush;
    const size_t cols_j_size =cols[j].size();
    for(size_t k=0; k<cols_j_size; k++){
      std::cerr << " " << cols[j][k] << std::flush;
    }
    std::cerr << std::endl;
  }
  
  //////// Elements ////////
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    //// m bound ////
    if(it->i>=m){
      ERROR_REPORT;
      std::cerr << "The indx of the " << it-data.begin() << "th element is out of the matrix size." << std::endl
                << "Its i index was " << it->i << "." << std::endl;
      exit(1);
    }
    
    //// n bound ////
    if(it->j>=n){
      ERROR_REPORT;
      std::cerr << "The indx of the " << it-data.begin() << "th element is out of the matrix size." << std::endl
                << "Its j index was " << it->j << "." << std::endl;
      exit(1);
    }
    
    //// double-listed ////
    for(std::vector<dcomponent>::const_iterator IT=it+1; IT!=data_end; IT++){
      if( it->i==IT->i && it->j==IT->j ){
        ERROR_REPORT;
        std::cerr << "The (" << it->i << ", " << it->j << ") component is double-listed at the " << it-data.begin() << "th and the" << IT-data.begin() << "the elements."<< std::endl;
        exit(1);
      }
    }
  }
  
  //////// ijc consistence ////////
  
  
  //////// NOTE ////////
  std::cerr << "# [NOTE]@dgsmatrix::checkup(): This sparse matrix is fine." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! swap two matrices */
inline void swap(dgsmatrix& A, dgsmatrix& B)
{CPPL_VERBOSE_REPORT;
  std::swap(A.n,B.n);
  std::swap(A.m,B.m);
  std::swap(A.data,B.data);
  std::swap(A.rows,B.rows);
  std::swap(A.cols,B.cols);
}

//=============================================================================
/*! convert user object to smart-temporary object */
inline _dgsmatrix _(dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  _dgsmatrix newmat;
  
  //////// shallow copy ////////
  newmat.n =mat.n;
  newmat.m =mat.m;
  std::swap(newmat.data,mat.data);
  std::swap(newmat.rows,mat.rows);
  std::swap(newmat.cols,mat.cols);
  
  //////// nullify ////////
  mat.m =0;
  mat.n =0;
  
  return newmat;
}
