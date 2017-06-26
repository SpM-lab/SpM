//=============================================================================
/*! operator() for non-const object */
inline double& dsymatrix::operator()(const CPPL_INT& i, const CPPL_INT& j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if( i >= j ) {
    //return array[i+n*j];
    return darray[j][i];
  } else {
    //return array[j+n*i];
    return darray[i][j];
  }
}

//=============================================================================
/*! operator() for const object */
inline double dsymatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if( i >= j ) {
    //return array[i+n*j];
    return darray[j][i];
  } else {
    //return array[j+n*i];
    return darray[i][j];
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! set value for const object */
inline dsymatrix& dsymatrix::set(const CPPL_INT& i, const CPPL_INT& j, const double& v) //const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if( i >= j ) {
    //array[i+n*j] = v;
    darray[j][i] =v;
  } else {
    //array[j+n*i] = v;
    darray[i][j] =v;
  }
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      s << " " << mat(i,j) << " ";
    }
    for(CPPL_INT j=i+1; j<mat.n; j++){
      s << "{" << mat(i,j) << "}";
    }
    s << std::endl;
  }
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void dsymatrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#dsymatrix " << n << std::endl;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      ofs << operator()(i,j) << " ";
    }
    ofs << std::endl;
  }
  
  ofs.close();
}

//=============================================================================
inline void dsymatrix::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s(filename);
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }

  std::string id;
  s >> id;
  if( id != "dsymatrix" && id != "#dsymatrix" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not dsymatrix." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  s >> n;
  resize(n);
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++ ){
      s >> operator()(i,j);
    }
  }
  if(s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << " ." << std::endl
              << "Most likely, there is not enough data components, or a linefeed code or space code is missing at the end of the last line." << std::endl;
    exit(1);
  }
  
  s >> id;
  if(!s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << "\"." << std::endl
              << "Most likely, there are extra data components." << std::endl;
    exit(1);
  }

  s.close();
}
