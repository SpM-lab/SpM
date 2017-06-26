//=============================================================================
/*! zcovector_small constructor */
template<CPPL_INT l>
inline zcovector_small<l>::zcovector_small()
{CPPL_VERBOSE_REPORT;
  ;
}

//=============================================================================
/*! zcovector_small constructor */
template<CPPL_INT l>
inline zcovector_small<l>::zcovector_small(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l!=vec.l ){
    ERROR_REPORT;
    std::cerr << "Vector sizes must be the same." << std::endl
              << "Your input was " << l << " and " << vec.l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT k=0; k<l; k++){
    array[k] =vec.array[k];
  }
}

//=============================================================================
/*! zcovector_small constructor */
template<CPPL_INT l>
inline zcovector_small<l>::zcovector_small(const comple& x, const comple& y)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l!=2 ){
    ERROR_REPORT;
    std::cerr << "The vector size must be 2." << std::endl
              << "The vector size you set was " << l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  array[0] =x;
  array[1] =y;
}

//=============================================================================
/*! zcovector_small constructor */
template<CPPL_INT l>
inline zcovector_small<l>::zcovector_small(const comple& x, const comple& y, const comple& z)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l!=3 ){
    ERROR_REPORT;
    std::cerr << "The vector size must be 3." << std::endl
              << "The vector size you set was " << l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  array[0] =x;
  array[1] =y;
  array[2] =z;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zcovector_small destructor */
template<CPPL_INT l>
inline zcovector_small<l>::~zcovector_small()
{CPPL_VERBOSE_REPORT;
  ;
}
