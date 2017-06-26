//=============================================================================
/*! operator() for object */
inline double& _dsymatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if( i >= j ){
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
inline std::ostream& operator<<(std::ostream& s, const _dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.n; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      if( i >= j ){
        s << " " << mat(i,j) << "  ";
      } else {
        s << "{" << mat(i,j) << "} ";
      }
    }
    s << std::endl;
  }
  
  mat.destroy();
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void _dsymatrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#dsymatrix " << n << std::endl;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++ ){
      ofs << operator()(i,j) << " ";
    }
    ofs << std::endl;
  }
  
  ofs.close();
  destroy();
}
