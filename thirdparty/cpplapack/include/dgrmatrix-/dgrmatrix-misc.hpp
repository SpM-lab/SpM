//=============================================================================
/*! clear all the matrix data and set the sizes 0 */
inline void dgrmatrix::clear()
{CPPL_VERBOSE_REPORT;
  m =0;
  n =0;
  a.clear();
  ia.clear();
  ja.clear();
}

//=============================================================================
/*! change the matrix into a zero matrix with no change in structure */
inline dgrmatrix& dgrmatrix::zero()
{CPPL_VERBOSE_REPORT;
  const std::vector<double>::iterator a_end =a.end();
  for(std::vector<double>::iterator it=a.begin(); it!=a_end; it++){
    (*it) =0.;
  }
  return *this;
}

//=============================================================================
/*! make a deep copy of the matrix */
inline void dgrmatrix::copy(const dgrmatrix& mat)
{CPPL_VERBOSE_REPORT;
  m =mat.m;
  n =mat.n;
  a =mat.a;
  ia =mat.ia;
  ja =mat.ja;
}

//=============================================================================
/*! check if the component is listed */
inline bool dgrmatrix::isListed(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  //////// search ////////
  int k_beg =ia[i]-1;
  int k_end =ia[i+1]-1;
  for(int k=k_beg; k<k_end; k++){
    if(ja[k]==j+1){
      return true;//found
    }
  }
  //////// not found ////////
  return false;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! health checkup */
inline void dgrmatrix::checkup()
{CPPL_VERBOSE_REPORT;
  //////// size check ////////
  if(m<0){
    ERROR_REPORT;
    std::cerr << "m<0" << std::endl;
    exit(1);    
  }
  if(n<0){
    ERROR_REPORT;
    std::cerr << "n<0" << std::endl;
    exit(1);    
  }
  if(a.size()!=ja.size()){
    ERROR_REPORT;
    std::cerr << "a.size()!=ja.size()" << std::endl;
    exit(1);
  }
  if(ia.size()!=size_t(m+1)){
    ERROR_REPORT;
    std::cerr << "ia.size()!=m+1" << std::endl;
    exit(1);
  }
  if(a.size()>size_t(m*n)){
    ERROR_REPORT;
    std::cerr << "a.size()>m*n" << std::endl;
    exit(1);
  }
  
  //////// index check ////////
  //not yet....
  
  std::cerr << "# [NOTE]@dgrmatrix::checkup(): This matrix is fine." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! swap two matrices */
inline void swap(dgrmatrix& A, dgrmatrix& B)
{CPPL_VERBOSE_REPORT;
  std::swap(A.n,B.n);
  std::swap(A.m,B.m);
  std::swap(A.a,B.a);
  std::swap(A.ia,B.ia);
  std::swap(A.ja,B.ja);
}
