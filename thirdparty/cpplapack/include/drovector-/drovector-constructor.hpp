//=============================================================================
/*! drovector constructor */
inline drovector::drovector()
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =0;
  cap =0;
  array =NULL;
}

//=============================================================================
/*! drovector copy constructor */
inline drovector::drovector(const drovector& vec)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =vec.l;
  cap =vec.cap;
  array =new double[cap];
  
  //////// copy ////////
  CPPL_INT inc =1;
  dcopy_(&l, vec.array, &inc, array, &inc);
}

//=============================================================================
/*! drovector constructor to cast _drovector */
inline drovector::drovector(const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =vec.l;
  cap =vec.cap;
  array =vec.array;
  
  vec.nullify();
}

//=============================================================================
/*! drovector constructor with size specification */
inline drovector::drovector(const CPPL_INT& _l, const CPPL_INT margin)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _l<0 || margin<0 ){
    ERROR_REPORT;
    std::cerr << "Vector size must be positive integers. " << std::endl
              << "Your input was (" << _l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// initialize ////////
  l =_l;
  cap =l+margin;
  array =new double[cap];
}

//=============================================================================
/*! drovector constructor with filename */
inline drovector::drovector(const char* filename)
{CPPL_VERBOSE_REPORT;
  array =NULL;
  read(filename);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! drovector destructor */
inline drovector::~drovector()
{CPPL_VERBOSE_REPORT;
  delete [] array;
}
