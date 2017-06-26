//=============================================================================
/*! _zhematrix constructor without arguments */
inline _zhematrix::_zhematrix()
  :m(n)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  n =0;
  array =NULL;
  darray =NULL;
}

//=============================================================================
/*! _zhematrix copy constructor */
inline _zhematrix::_zhematrix(const _zhematrix& mat)
  :m(n)
{CPPL_VERBOSE_REPORT;
  //////// initialize ////////
  n =mat.n;
  array =mat.array;
  darray =mat.darray;
  
  mat.nullify();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix destructor */
inline _zhematrix::~_zhematrix()
{CPPL_VERBOSE_REPORT;
  delete [] array;
  delete [] darray;
}
