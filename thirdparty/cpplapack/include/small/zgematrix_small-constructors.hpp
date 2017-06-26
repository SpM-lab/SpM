//=============================================================================
/*! zgematrix_small constructor */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>::zgematrix_small()
{CPPL_VERBOSE_REPORT;
  ;
}

//=============================================================================
/*! zgematrix_small constructor */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>::zgematrix_small(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( m!=mat.m || n!=mat.n ){
    ERROR_REPORT;
    std::cerr << "Matrix sizes must be the same." << std::endl
              << "Your input was " << m << "x" << n << " and " << mat.m << "x" << mat.n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT k=0; k<m*n; k++){
    array[k] =mat.array[k];
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix_small destructor */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>::~zgematrix_small()
{CPPL_VERBOSE_REPORT;
  ;
}
