//=============================================================================
/*! dcovector constructor */
inline dcovector::dcovector()
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =0;
  cap =0;
  array =NULL;
}

//=============================================================================
/*! dcovector copy constructor */
inline dcovector::dcovector(const dcovector& vec)
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
/*! dcovector constructor to cast _dcovector */
inline dcovector::dcovector(const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  l =vec.l;
  cap =vec.cap;
  array =vec.array;
  
  vec.nullify();
}

//=============================================================================
/*! dcovector constructor with size specification */
inline dcovector::dcovector(const CPPL_INT& _l, const CPPL_INT margin)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _l<0 ){
    ERROR_REPORT;
    std::cerr << "Vector size must be positive integers. " << std::endl
              << "Your input was (" << _l << ")." << std::endl;
    exit(1);
  }
  if( margin<0 ){
    ERROR_REPORT;
    std::cerr << "Vector margin must be zero or above. " << std::endl
              << "Your input was (" << _l << ", " << margin << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// initialize ////////
  l =_l;
  cap =l+margin;
  array =new double[cap];
}

//=============================================================================
/*! dcovector constructor with filename */
inline dcovector::dcovector(const char* filename)
{CPPL_VERBOSE_REPORT;
  array =NULL;
  read(filename);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector destructor */
inline dcovector::~dcovector()
{CPPL_VERBOSE_REPORT;
  //////// delete array ////////
  delete [] array;
}
