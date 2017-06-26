//=============================================================================
/*! zhematrix constructor without arguments */
inline zhematrix::zhematrix()
  : m(n)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  n = 0;
  array =NULL;
  darray =NULL;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix copy constructor */
inline zhematrix::zhematrix(const zhematrix& mat)
  : m(n)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  n =mat.n;
  array =new comple[n*n];
  darray =new comple*[n];
  for(int i=0; i<n; i++){
    darray[i] =&array[i*n];
  }
  
  //////// copy ////////
  CPPL_INT size =n*n;
  CPPL_INT inc =1;
  zcopy_(&size, mat.array, &inc, array, &inc);
}

//=============================================================================
/*! zhematrix constructor to cast _zhematrix */
inline zhematrix::zhematrix(const _zhematrix& mat)
  : m(n)
{CPPL_VERBOSE_REPORT;
  n =mat.n;
  array =mat.array;
  darray =mat.darray;
  
  mat.nullify();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix constructor with size specification */
inline zhematrix::zhematrix(const CPPL_INT& _n)
  : m(n)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _n<0 ){
    ERROR_REPORT;
    std::cerr << "Matrix sizes must be positive integers. " << std::endl
              << "Your input was (" << _n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// initialize ////////
  n =_n;
  array =new comple[n*n];
  darray =new comple*[n];
  for(int i=0; i<n; i++){
    darray[i] =&array[i*n];
  }
}

//=============================================================================
/*! zhematrix constructor with filename */
inline zhematrix::zhematrix(const char* filename)
  : m(n)
{CPPL_VERBOSE_REPORT;
  array =NULL;
  darray =NULL;
  
  //// read ////
  read(filename);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix destructor */
inline zhematrix::~zhematrix()
{CPPL_VERBOSE_REPORT;
  //////// delete array ////////
  delete [] array;
  delete [] darray;
}
