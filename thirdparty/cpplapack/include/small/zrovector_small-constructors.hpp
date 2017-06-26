//=============================================================================
/*! zrovector_small constructor */
template<CPPL_INT l>
inline zrovector_small<l>::zrovector_small()
{CPPL_VERBOSE_REPORT;
  ;
}

//=============================================================================
/*! zrovector_small constructor */
template<CPPL_INT l>
inline zrovector_small<l>::zrovector_small(const zrovector& vec)
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
/*! zrovector_small constructor */
template<CPPL_INT l>
inline zrovector_small<l>::zrovector_small(const comple& x, const comple& y)
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
/*! zrovector_small constructor */
template<CPPL_INT l>
inline zrovector_small<l>::zrovector_small(const comple& x, const comple& y, const comple& z)
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
/*! zrovector_small destructor */
template<CPPL_INT l>
inline zrovector_small<l>::~zrovector_small()
{CPPL_VERBOSE_REPORT;
  ;
}
