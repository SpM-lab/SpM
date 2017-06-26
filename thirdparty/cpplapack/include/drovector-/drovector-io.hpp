//=============================================================================
/*! operator() for non-const object */
inline double& drovector::operator()(const CPPL_INT& i)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || l<=i ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the vector size." << std::endl
              << "Your input is (" << i << "), whereas the vector size is " << l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  return array[i];
}

//=============================================================================
/*! operator() for const object */
inline double drovector::operator()(const CPPL_INT& i) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || l<=i ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the vector size." << std::endl
              << "Your input is (" << i << "), whereas the vector size is " << l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  return array[i];
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! set value for const object */
inline drovector& drovector::set(const CPPL_INT& i, const double& v) //const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || l<=i ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the vector size." << std::endl
              << "Your input is (" << i << "), whereas the vector size is " << l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  array[i] =v;
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const drovector& vec)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<vec.l; i++){ s << " " << vec.array[i]; }
  s << std::endl;
  
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void drovector::write(const char *filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#drovector" << " " << l << std::endl;
  for(CPPL_INT i=0; i<l; i++){
    ofs << operator()(i) << " ";
  }
  ofs << std::endl;
  
  ofs.close();
}

//=============================================================================
inline void drovector::read(const char *filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }

  std::string id;
  s >> id;
  if( id != "drovector" && id != "#drovector" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not drovector." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  s >> l;
  resize(l);
  for(CPPL_INT i=0; i<l; i++) { s >> operator()(i); }
  if(s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << "\"." << std::endl
              << "Most likely, there is a lack of data components, or a linefeed code or space code is missing at the end of the last line." << std::endl;
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
