//=============================================================================
/*! dsymatrix_small constructor */
template<CPPL_INT n>
inline dsymatrix_small<n>::dsymatrix_small()
{CPPL_VERBOSE_REPORT;
  ;
}

//=============================================================================
/*! dsymatrix_small constructor */
template<CPPL_INT n>
inline dsymatrix_small<n>::dsymatrix_small(const dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( n!=mat.n ){
    ERROR_REPORT;
    std::cerr << "Matrix sizes must be the same." << std::endl
              << "Your input was " << n << "x" << n << " and " << mat.m << "x" << mat.n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT k=0; k<(n*(n+1))/2; k++){
    array[k] =mat.array[k];
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dsymatrix_small destructor */
template<CPPL_INT n>
inline dsymatrix_small<n>::~dsymatrix_small()
{CPPL_VERBOSE_REPORT;
  ;
}
