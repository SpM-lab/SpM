//=============================================================================
/*! clear all the matrix data and set the sizes 0 */
inline void dssmatrix::clear()
{CPPL_VERBOSE_REPORT;
  n =0;
  data.clear();
  line.clear();
}

//=============================================================================
/*! change the matrix into a zero matrix */
inline dssmatrix& dssmatrix::zero()
{CPPL_VERBOSE_REPORT;
  data.resize(0);
  for(CPPL_INT i=0; i<n; i++){ line[i].resize(0); }
  return *this;
}

//=============================================================================
/*! change the matrix into an identity matrix */
inline dssmatrix& dssmatrix::identity()
{CPPL_VERBOSE_REPORT;
  zero();
  for(CPPL_INT i=0; i<n; i++){
    put(i,i,1.);
  }
  return *this;
}

//=============================================================================
/*! change sign(+/-) of the matrix */
inline void dssmatrix::chsign()
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator data_end =data.end();
  for(std::vector<dcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v =-it->v;
  }
}

//=============================================================================
/*! make a deep copy of the matrix */
inline void dssmatrix::copy(const dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
  n =mat.n;
  data =mat.data;
  line =mat.line;
}

//=============================================================================
/*! make a shallow copy of the matrix\n
  This function is not designed to be used in project codes. */
inline void dssmatrix::shallow_copy(const _dssmatrix& mat)
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
inline dssmatrix& dssmatrix::resize(const CPPL_INT& _n, const CPPL_INT _c, const CPPL_INT _l)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _n<0 || _c<0 || _l<0 ){
    ERROR_REPORT;
    std::cerr << "Matrix sizes, the length of arrays, and line size must be positive integers. " << std::endl
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
inline void dssmatrix::stretch(const CPPL_INT& dn)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( n+dn<0 ){
    ERROR_REPORT;
    std::cerr << "The new matrix size must be larger than zero." << std::endl
              << "Your input was (" << dn << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// zero ////////
  if(dn==0){ return; }

  //////// non-zero ////////
  n +=dn;
  
  if(dn<0){
    //// delete components over the new size ////
    const std::vector<dcomponent>::reverse_iterator data_rend =data.rend();
    for(std::vector<dcomponent>::reverse_iterator it=data.rbegin(); it!=data_rend; it++){
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
inline bool dssmatrix::isListed(const CPPL_INT& i, const CPPL_INT& j) const
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
    if(data[*p].i==ii && data[*p].j==jj){ return 1; }
  }
  
  return 0;
}

//=============================================================================
/*! return the element number of the component */
inline CPPL_INT dssmatrix::number(const CPPL_INT& i, const CPPL_INT& j) const
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
    if(data[*p].i==ii && data[*p].j==jj){ return *p; }
  }
  
  return -1;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! get row of the matrix */
inline _drovector dssmatrix::row(const CPPL_INT& _m) const
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
  
  const std::vector<CPPL_INT>::const_iterator line__m_end =line[_m].end();
  for(std::vector<CPPL_INT>::const_iterator p=line[_m].begin(); p!=line__m_end; p++){
    if(data[*p].i==_m){
      vec(data[*p].j) =data[*p].v;
    }
    else{
      vec(data[*p].i) =data[*p].v;
    }
  }
  return _(vec);
}

//=============================================================================
/*! get column of the matrix */
inline _dcovector dssmatrix::col(const CPPL_INT& _n) const
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
  
  const std::vector<CPPL_INT>::const_iterator line__n_end =line[_n].end();
  for(std::vector<CPPL_INT>::const_iterator p=line[_n].begin(); p!=line__n_end; p++){
    if(data[*p].i==_n){
      vec(data[*p].j) =data[*p].v;
    }
    else{
      vec(data[*p].i) =data[*p].v;
    }
  }
  return _(vec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! erase components less than DBL_MIN */
inline void dssmatrix::diet(const double eps)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::reverse_iterator data_rend =data.rend();
  for(std::vector<dcomponent>::reverse_iterator it=data.rbegin(); it!=data_rend; it++){
    if( fabs(it->v)<eps ){ del( CPPL_INT(data_rend-it-1) ); }
  }
}

//=============================================================================
/*! reorder components so that all diagonal componets are placed in front */
/*
inline CPPL_INT dssmatrix::diag_front()
{CPPL_VERBOSE_REPORT;
  //////// set initial dsize ////////
  CPPL_INT dsize(0);
  const std::vector<dcomponent>::iterator data_end =data.end();
  for(std::vector<dcomponent>::iterator it=data.begin(); it!=data.end(); it++){
    if(it->i==it->j){ dsize++; }
    else{ break; }
  }
  
  //////// swapping loop ////////
  for(std::vector<dcomponent>::reverse_iterator it=data.rbegin(); it!=data.rend()-dsize; it++){
    if(it->i==it->j){//is diag
      CPPL_INT c(data.rend()-it-1);//current it's index
      CPPL_INT i(data[dsize].i), j(data[dsize].j), k(it->i);
      //// search (k,k) line ////
      for(std::vector<CPPL_INT>::iterator p=line[k].begin(); p!=line[k].end(); p++){
        if(CPPL_INT(data[*p].i)==k && CPPL_INT(data[*p].j)==k){ *p=dsize; }
      }
      //// search (i,j) line ////
      for(std::vector<CPPL_INT>::iterator p=line[i].begin(); p!=line[i].end(); p++){
        if(CPPL_INT(*p)==dsize){ *p=c; }
      }
      //// search (j,i) line ////
      if(i!=j){
        for(std::vector<CPPL_INT>::iterator p=line[j].begin(); p!=line[j].end(); p++){
          if(CPPL_INT(*p)==dsize){ *p=c; }
        }
      }
      else{//i==j
        it--;
      }
      //// swap ////
      std::swap(data[dsize],data[c]);
      //// update ////
      dsize++;
    }
  }
  
  return dsize;
}
*/

//=============================================================================
/*! reorder components */
inline void dssmatrix::reorder(const bool mode)
{CPPL_VERBOSE_REPORT;
  //// sort data ////
  if(mode==0){
    std::sort(data.begin(), data.end(), dcomponent::ilt);
  }
  else{
    std::sort(data.begin(), data.end(), dcomponent::jlt);
  }
  //// rebuild line ////
  rebuild();
}

//=============================================================================
/*! rebuild line */
inline void dssmatrix::rebuild()
{CPPL_VERBOSE_REPORT;
  //// clear line ////
  for(CPPL_INT i=0; i<n; i++){ line[i].resize(0); }
  
  //// build line ////
  CPPL_INT c(0);
  const std::vector<dcomponent>::iterator data_end =data.end();
  for(std::vector<dcomponent>::iterator it=data.begin(); it!=data_end; it++){
    line[it->i].push_back(c);
    if( (it->i) != (it->j) ){
      line[it->j].push_back(c);
    }
    c++;
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! health checkup */
inline void dssmatrix::checkup()
{CPPL_VERBOSE_REPORT;
  //////// write ////////
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    std::cerr << "array[" << it-data.begin() << "] = (" << it->i << "," << it->j << ") = " << it->v << std::endl;
  }
  std::cerr << std::endl;
  
  for(CPPL_INT i=0; i<n; i++){
    std::cerr << "line[" << i << "] =" << std::flush;
    const size_t line_i_size =line[i].size();
    for(size_t k=0; k<line_i_size; k++){
      std::cerr << " " << line[i][k] << std::flush;
    }
    std::cerr << std::endl;
  }
  std::cerr << std::endl;
  
  //////// Elements ////////
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    //// m bound ////
    if(it->i>=n){
      ERROR_REPORT;
      std::cerr << "The indx of the " << it-data.begin() << "th element is out of the matrix size." << std::endl
                << "Its i index was " << it->i << "." << std::endl;
      exit(1);
    }
    
    //// n bound ////
    if(it->j>=n){
      ERROR_REPORT;
      std::cerr << "The jndx of the " << it-data.begin() << "th element is out of the matrix size." << std::endl
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
  
  //////// NOTE ////////
  std::cerr << "# [NOTE]@dssmatrix::checkup(): This symmetric sparse matrix is fine." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! swap two matrices */
inline void swap(dssmatrix& A, dssmatrix& B)
{CPPL_VERBOSE_REPORT;
  std::swap(A.n,B.n);
  std::swap(A.data,B.data);
  std::swap(A.line,B.line);
}

//=============================================================================
/*! convert user object to smart-temporary object */
inline _dssmatrix _(dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
  _dssmatrix newmat;
  
  //////// shallow copy ////////
  newmat.n =mat.n;
  std::swap(newmat.data,mat.data);
  std::swap(newmat.line,mat.line);
  
  //////// nullify ////////
  mat.n =0;
  
  return newmat;
}
