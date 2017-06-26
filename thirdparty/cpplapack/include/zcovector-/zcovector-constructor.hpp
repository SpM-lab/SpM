//=============================================================================
/*! zcovector constructor */
inline zcovector::zcovector()
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =0;
  array =NULL;
}

//=============================================================================
/*! zcovector copy constructor */
inline zcovector::zcovector(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =vec.l;
  array =new comple[l];
  
  //////// copy ////////
  CPPL_INT inc =1;
  zcopy_(&l, vec.array, &inc, array, &inc);
}

//=============================================================================
/*! zcovector constructor to cast _zcovector */
inline zcovector::zcovector(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =vec.l;
  array =vec.array;
  
  vec.nullify();
}

//=============================================================================
/*! zcovector constructor with size specification */
inline zcovector::zcovector(const CPPL_INT& _l)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _l<0 ){
    ERROR_REPORT;
    std::cerr << "Vector size must be positive integers. " << std::endl
              << "Your input was (" << _l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// initialize ////////
  l =_l;
  array =new comple[l];
}

//=============================================================================
/*! zcovector constructor with filename */
inline zcovector::zcovector(const char* filename)
{CPPL_VERBOSE_REPORT;
  array =NULL;
  
  //// copy ////
  read(filename);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zcovector destructor */
inline zcovector::~zcovector()
{CPPL_VERBOSE_REPORT;
  //////// delete array ////////
  delete [] array;
}
