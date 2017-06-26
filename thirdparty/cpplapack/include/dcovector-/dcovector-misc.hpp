//=============================================================================
/*! clear vector */
inline void dcovector::clear()
{CPPL_VERBOSE_REPORT;
  l =0;
  cap =0;
  delete [] array;
  array =NULL;
}

//=============================================================================
/*! make vector into zero vector */
inline dcovector& dcovector::zero()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){ array[i] =0.0; }
  return *this;
}

//=============================================================================
/*! change sign(+/-) of the vector */
inline void dcovector::chsign()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){ array[i] =-array[i]; }
}

//=============================================================================
/*! make a deep copy of the dcovector */
inline void dcovector::copy(const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  cap =vec.cap;
  delete [] array;
  array =new double[cap];
  CPPL_INT inc =1;
  dcopy_(&vec.l, vec.array, &inc, array, &inc);
}

//=============================================================================
/*! make a shallow copy of the vector\n
 This function is not desinged to be used in project codes. */
inline void dcovector::shallow_copy(const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  cap =vec.cap;
  delete [] array;
  array =vec.array;
  
  vec.nullify();
}

//=============================================================================
/*! make an alias of the vector\n
 Be carefull to use this function not to cause double free. */
inline void dcovector::alias(const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  cap =vec.cap;
  delete [] array;
  array =vec.array;
}

//=============================================================================
/*! unalias the vector */
inline void dcovector::unalias()
{CPPL_VERBOSE_REPORT;
  l =0;
  cap =0;
  array =NULL;
}

//=============================================================================
/*! resize vector */
inline dcovector& dcovector::resize(const CPPL_INT& _l, const CPPL_INT margin)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _l<0 ){
    ERROR_REPORT;
    std::cerr << "Vector size must be positive integers." << std::endl
              << "Your input was (" << _l << ", " << margin << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  l =_l;
  cap =l+margin;
  delete [] array;
  array =new double[cap];
  
  return *this;
}

//=============================================================================
/*! stretch or shrink vector */
inline void dcovector::stretch(const CPPL_INT& dl)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l+dl<0 ){
    ERROR_REPORT;
    std::cerr << "Vector size must be positive integers." << std::endl
              << "Your input was (" << dl << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// zero ////////
  if(dl==0){ return; }
  
  //////// non-zero ////////
  l +=dl;
  if(l>cap){
    while(l>cap){
      cap++;
      cap*=2;
    }
    CPPL_INT newl =l-dl;
    CPPL_INT inc =1;
    double* newArray(new double[cap]);
    dcopy_(&newl, array, &inc, newArray, &inc);
    delete [] array;
    array =newArray;
  }
}

//=============================================================================
/*! swap two vectors */
inline void swap(dcovector& u, dcovector& v)
{CPPL_VERBOSE_REPORT;
  CPPL_INT u_cap(u.cap), u_l(u.l);
  double* u_array(u.array);
  u.l=v.l; u.cap=v.cap; u.array=v.array;
  v.l=u_l; v.cap=u_cap; v.array=u_array;
}

//=============================================================================
/*! convert user object to smart-temporary object */
inline _dcovector _(dcovector& vec)
{CPPL_VERBOSE_REPORT;
  _dcovector newvec;
  
  //////// shallow copy ////////
  newvec.l =vec.l;
  newvec.cap =vec.cap;
  newvec.array =vec.array;
  
  //////// nullify ////////
  vec.l =0;
  vec.cap =0;
  vec.array =NULL;
  
  return newvec;
}
