//=============================================================================
/*! dsymatrix constructor without arguments */
inline dsymatrix::dsymatrix()
  : m(n)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  n = 0;
  array =NULL;
  darray =NULL;
}

//=============================================================================
/*! dsymatrix copy constructor */
inline dsymatrix::dsymatrix(const dsymatrix& mat)
  : m(n)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  n =mat.n;
  array =new double[n*n];
  darray =new double*[n];
  for(int i=0; i<n; i++){
    darray[i] =&array[i*n];
  }
  
  //////// copy ////////
  CPPL_INT size =n*n;
  CPPL_INT inc =1;
  dcopy_(&size, mat.array, &inc, array, &inc);
}

//=============================================================================
/*! dsymatrix constructor to cast _dsymatrix */
inline dsymatrix::dsymatrix(const _dsymatrix& mat)
  : m(n)
{CPPL_VERBOSE_REPORT;
  n =mat.n;
  array =mat.array;
  darray =mat.darray;
  
  mat.nullify();
}

//=============================================================================
/*! dsymatrix copy constructor to cast dssmatrix */
/*
inline dsymatrix::dsymatrix(const dssmatrix& mat)
  : m(n)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  n =mat.n;
  array =new double[n*n];
  darray =new double*[n];
  for(int i=0; i<n; i++){ darray[i] =&array[i*n]; }
  
  //////// copy ////////
  zero();
  for(int c=0; c<mat.vol; c++){
    (*this)(mat.indx[c],mat.jndx[c]) =mat.array[c];
  }
}
*/

//=============================================================================
/*! dsymatrix constructor to cast _dssmatrix */
/*
inline dsymatrix::dsymatrix(const _dssmatrix& mat)
  : m(n)
{CPPL_VERBOSE_REPORT;
  n =mat.n;
  array =new double[n*n];
  darray =new double*[n];
  for(int i=0; i<n; i++){ darray[i] =&array[i*n]; }
  
  //////// copy ////////
  zero();
  for(int c=0; c<mat.vol; c++){
    (*this)(mat.indx[c],mat.jndx[c]) =mat.array[c];
  }
  
  mat.nullify();
}
*/


//=============================================================================
/*! dsymatrix constructor with size specification */
inline dsymatrix::dsymatrix(const CPPL_INT& _n)
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
  array =new double[n*n];
  darray =new double*[n];
  for(int i=0; i<n; i++){
    darray[i] =&array[i*n];
  }
}

//=============================================================================
/*! dsymatrix constructor with filename */
inline dsymatrix::dsymatrix(const char* filename)
  : m(n)
{CPPL_VERBOSE_REPORT;
  array =NULL;
  darray =NULL;
  
  //// copy ////
  read(filename);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dsymatrix destructor */
inline dsymatrix::~dsymatrix()
{CPPL_VERBOSE_REPORT;
  delete [] array;
  delete [] darray;
}
