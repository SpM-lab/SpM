//=============================================================================
/*! operator() for const object */
inline double& _dgbmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j || i-j>kl || j-i>ku ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is (" << m << "," << n << " , " << kl << "," << ku << ")." << std::endl;
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
inline std::ostream& operator<<(std::ostream& s, const _dgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      if( i-j>mat.kl || j-i>mat.ku ){ s << " x"; }
      else{ s << " " << mat(i,j); }
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
inline void _dgbmatrix::write(const char *filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#dgbmatrix" << " " << m << " " << n << " " << kl << " " << ku << std::endl;
  for(CPPL_INT i=0; i<m; i++){
    const CPPL_INT jmax =std::min(n,i+ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-kl); j<jmax; j++){
      ofs << operator()(i,j) << " ";
    }
    ofs << std::endl;
  }
  
  ofs.close();
  destroy();
}
