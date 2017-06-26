//=============================================================================
/*! zrovector constructor */
inline zrovector::zrovector()
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =0;
  array =NULL;
}

//=============================================================================
/*! zrovector copy constructor */
inline zrovector::zrovector(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =vec.l;
  array =new comple[l];
  
  //////// copy ////////
  CPPL_INT inc =1;
  zcopy_(&l, vec.array, &inc, array, &inc);
}

//=============================================================================
/*! zrovector constructor to cast _zrovector */
inline zrovector::zrovector(const _zrovector& vec)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =vec.l;
  array =vec.array;
  
  vec.nullify();
}

//=============================================================================
/*! zrovector constructor with size specification */
inline zrovector::zrovector(const CPPL_INT& _l)
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
/*! zrovector constructor with filename */
inline zrovector::zrovector(const char* filename)
{CPPL_VERBOSE_REPORT;
  array =NULL;
  
  //// read ////
  read(filename);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector destructor */
inline zrovector::~zrovector()
{CPPL_VERBOSE_REPORT;
  //////// delete array ////////
  delete [] array;
}
