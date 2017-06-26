//=============================================================================
/*! complete the upper-right components */
inline void dsymatrix::complete() const
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<i; j++){
      darray[i][j] =darray[j][i];
    }
  }
}

//=============================================================================
/*! clear all the matrix data and set the sizes 0 */
inline void dsymatrix::clear()
{CPPL_VERBOSE_REPORT;
  n =0;
  delete [] array;
  array =NULL;
  delete [] darray;
  darray =NULL;
}

//=============================================================================
/*! change the matrix into a zero matrix */
inline dsymatrix& dsymatrix::zero()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      darray[j][i] =0.;
    }
  }
  return *this;
}

//=============================================================================
/*! change the matrix into an identity matrix */
inline dsymatrix& dsymatrix::identity()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<n; j++){
    darray[j][j] =1.;
    for(CPPL_INT i=j+1; i<n; i++){
      darray[j][i] =0.;
    }
  }
  return *this;
}

//=============================================================================
/*! change sign(+/-) of the matrix */
inline void dsymatrix::chsign()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      darray[j][i] =-darray[j][i];
    }
  }
}

//=============================================================================
/*! make a deep copy of the matrix */
inline void dsymatrix::copy(const dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  n =mat.n;
  delete [] array;
  array =new double[n*n];
  delete [] darray;
  darray =new double*[n];
  for(int i=0; i<n; i++){
    darray[i] =&array[i*n];
  }
  
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=j; i<n; i++){
      darray[j][i] =mat.darray[j][i];
    }
  }
}

//=============================================================================
/*! make a shallow copy of the matrix\n
  This function is not designed to be used in project codes. */
inline void dsymatrix::shallow_copy(const _dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  n =mat.n;
  delete [] array;
  array =mat.array;
  delete [] darray;
  darray =mat.darray;
  
  mat.nullify();
}

//=============================================================================
/*! resize the matrix */
inline dsymatrix& dsymatrix::resize(const CPPL_INT& _n)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _n<0 ){
    ERROR_REPORT;
    std::cerr << "Matrix sizes must be positive integers." << std::endl
              << "Your input was (" << _n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  n =_n;
  delete [] array;
  array =new double[n*n];
  delete [] darray;
  darray =new double*[n];
  for(int i=0; i<n; i++){
    darray[i] =&array[i*n];
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! get row of the matrix */
inline _drovector dsymatrix::row(const CPPL_INT& _m) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _m<0 || _m>m ){
    ERROR_REPORT;
    std::cerr << "Input row number must be between 0 and " << m << "." << std::endl
              << "Your input was " << _m << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  drovector v(n);
  for(CPPL_INT j=0; j<n; j++){
    v(j)=(*this)(_m,j);
  }
  return _(v);
}

//=============================================================================
/*! get column of the matrix */
inline _dcovector dsymatrix::col(const CPPL_INT& _n) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _n<0 || _n>n ){
    ERROR_REPORT;
    std::cerr << "Input row number must be between 0 and " << n << "." << std::endl
              << "Your input was " << _n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dcovector v(m);
  for(CPPL_INT i=0; i<m; i++){
    v(i)=(*this)(i,_n);
  }
  return _(v);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! swap two matrices */
inline void swap(dsymatrix& A, dsymatrix& B)
{CPPL_VERBOSE_REPORT;
  CPPL_INT A_n =A.n;
  double* A_array =A.array;
  double** A_darray =A.darray;
  A.n=B.n; A.array=B.array; A.darray=B.darray;
  B.n=A_n; B.array=A_array; B.darray=A_darray;
}

//=============================================================================
/*! convert user object to smart-temporary object */
inline _dsymatrix _(dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  _dsymatrix newmat;

  //////// shallow copy ////////
  newmat.n =mat.n;
  newmat.array =mat.array;
  newmat.darray =mat.darray;
  
  //////// nullify ////////
  mat.n =0;
  mat.array =NULL;
  mat.darray =NULL;
  
  return newmat;
}
