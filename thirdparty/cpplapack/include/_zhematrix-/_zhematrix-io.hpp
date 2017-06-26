//=============================================================================
/*! operator() for object */
inline zhecomplex _zhematrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is (" << n << "," << n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if(i>=j){ return zhecomplex(i,j, darray[j][i]); }
  else    { return zhecomplex(i,j, darray[i][j]); }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const _zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.n; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      if(i>j){
        s << " " << mat(i,j) << " ";
      }
      else if(i==j){
        s << " " << std::real(mat(i,i)) << " ";
      }
      else{
        s << "{" << std::conj(mat.darray[i][j]) << "} ";
      }
    }
    s << std::endl;
    
#ifdef  CPPL_DEBUG
    if(std::fabs(std::imag(mat(i,i))) > DBL_MIN){
      WARNING_REPORT;
      std::cerr << "The " << i << "th diagonal component of the zhematrix is not a real number." << std::endl;
    }
#endif//CPPL_DEBUG
  }
  
  mat.destroy();
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void _zhematrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#zhematrix " << n << std::endl;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++ ){
      ofs << operator()(i,j) << " ";
    }
    ofs << std::endl;
    
#ifdef  CPPL_DEBUG
    if(std::fabs(std::imag(operator()(i,i))) > DBL_MIN){
      WARNING_REPORT;
      std::cerr << "The " << i << "th diagonal component of the zhematrix is not a real number." << std::endl;
    }
#endif//CPPL_DEBUG
  }
  
  ofs.close();
  destroy();
}
