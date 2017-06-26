//=============================================================================
/*! operator() for non-const object */
inline comple& zgbmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j || i-j>kl || j-i>ku ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << " with kl=" << kl << ", ku=" << ku << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //return array[ku+i+(kl+ku)*j];
  return darray[j][ku-j+i];
}

//=============================================================================
/*! operator() for const object */
inline comple zgbmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j || i-j>kl || j-i>ku ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << " with kl=" << kl << ", ku=" << ku << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //return array[ku+i+(kl+ku)*j];
  return darray[j][ku-j+i];
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! set value for const object */
inline zgbmatrix& zgbmatrix::set(const CPPL_INT& i, const CPPL_INT& j, const comple& v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j || i-j>kl || j-i>ku ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << " with kl=" << kl << ", ku=" << ku << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //array[ku+i+(kl+ku)*j] =v;
  darray[j][ku-j+i] =v;
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      if( i-j>mat.kl || j-i>mat.ku ){ s << " x"; }
      else{ s << " " << mat(i,j); }
    }
    s << std::endl;
  }
  
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void zgbmatrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#zgbmatrix" << " " << m << " " << n << " " << kl << " " << ku << std::endl;
  for(CPPL_INT i=0; i<m; i++){
    const CPPL_INT jmax =std::min(n,i+ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-kl); j<jmax; j++){
      ofs << operator()(i,j) << " ";
    }
    ofs << std::endl;
  }
  
  ofs.close();
}

//=============================================================================
inline void zgbmatrix::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }

  std::string id;
  s >> id;
  if( id != "zgbmatrix" && id != "#zgbmatrix" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not zgbmatrix." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  s >> m >> n >> kl >> ku;
  resize(m, n, kl, ku);
  for(CPPL_INT i=0; i<m; i++){
    const CPPL_INT jmax =std::min(n,i+ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-kl); j<jmax; j++){
      s >> operator()(i,j);
    }
  }
  if(s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << "\"." << std::endl
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
