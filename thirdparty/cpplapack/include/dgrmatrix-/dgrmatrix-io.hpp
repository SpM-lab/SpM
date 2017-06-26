//=============================================================================
/*! operator() for const object */
inline double dgrmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// search (i,j) component ////
  int k_beg =ia[i]-1;
  int k_end =ia[i+1]-1;
  for(int k=k_beg; k<k_end; k++){
    if(j==ja[k]-1){
      return a[k];
    }
  }
  
  //// (i,j) component was not found ////
  return 0.0;
}

//=============================================================================
/*! operator() for const object */
inline double& dgrmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// search (i,j) component ////
  int k_beg =ia[i]-1;
  int k_end =ia[i+1]-1;
  for(int k=k_beg; k<k_end; k++){
    if(j==ja[k]-1){
      return a[k];
    }
  }
  
  //// (i,j) component was not found ////
  ERROR_REPORT;
  std::cerr << "dgrmatrix does not allow component addition with operator()." << std::endl;
  exit(1);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const dgrmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m; i++){
    int k_beg =mat.ia[i]-1;
    int k_end =mat.ia[i+1]-1;
    int j =0;
    for(int k=k_beg; k<k_end; k++){
      if(j<mat.ja[k]-1){
        for(; j<mat.ja[k]-1; j++){
          s << " x";
        }
      }
      s << " " << mat.a[k];
      j++;
    }
    for(; j<mat.n; j++){
      s << " x";
    }    
    s << std::endl;
  }
  
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void dgrmatrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#dgrmatrix" << " " << m << " " << n << " " << a.size() << std::endl;
  
  size_t a_size =a.size();
  for(size_t k=0; k<a_size; k++){
    ofs << " " << a[k];
  }
  ofs << "\n";
  
  size_t ia_size =ia.size();
  for(size_t k=0; k<ia_size; k++){
    ofs << " " << ia[k];
  }
  ofs << "\n";
  
  size_t ja_size =ja.size();
  for(size_t k=0; k<ja_size; k++){
    ofs << " " << ja[k];
  }
  ofs << "\n" << std::flush;
  
  ofs.close();
}

//=============================================================================
inline void dgrmatrix::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }
  
  std::string id;
  s >> id;
  if( id != "dgrmatrix" && id != "#dgrmatrix" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not dgrmatrix." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  //////// read ////////
  size_t a_size;
  s >> m >> n >> a_size;
  a.resize(a_size);
  ia.resize(m+1);
  ja.resize(a_size);
  
  for(size_t k=0; k<a_size; k++){
    s >> a[k];
    if(s.fail()){
      ERROR_REPORT;
      exit(1);
    }
  }
  for(CPPL_INT k=0; k<=m; k++){
    s >> ia[k];
    if(s.fail()){
      ERROR_REPORT;
      exit(1);
    }
  }
  for(size_t k=0; k<a_size; k++){
    s >> ja[k];
    if(s.fail()){
      ERROR_REPORT;
      exit(1);
    }
  }
  
  //////// garbage ////////
  s >> a_size;
  if(!s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << " ." << std::endl;
    exit(1);
  }
  s.close();
}
