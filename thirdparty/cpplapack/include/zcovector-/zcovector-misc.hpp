//=============================================================================
/*! clear vector */
inline void zcovector::clear()
{CPPL_VERBOSE_REPORT;
  l =0;
  delete [] array;
  array =NULL;
}

//=============================================================================
/*! make vector into zero vector */
inline zcovector& zcovector::zero()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    array[i] =comple(0.,0.);
  }
  return *this;
}

//=============================================================================
/*! change sign(+/-) of the vector */
inline void zcovector::chsign()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    array[i] =-array[i];
  }
}

//=============================================================================
/*! make a deep copy of the zcovector */
inline void zcovector::copy(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  delete [] array;
  array =new comple[vec.l];
  
  CPPL_INT inc =1;
  zcopy_(&vec.l, vec.array, &inc, array, &inc);
}

//=============================================================================
/*! make a shallow copy of the vector\n
 This function is not desinged to be used in project codes. */
inline void zcovector::shallow_copy(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  delete [] array;
  array =vec.array;
  
  vec.nullify();
}

//=============================================================================
/*! make an alias of the vector\n
  Be carefull to use this function not to cause double free. */
inline void zcovector::alias(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  //cap =vec.cap;
  delete [] array;
  array =vec.array;
}

//=============================================================================
/*! unalias the vector */
inline void zcovector::unalias()
{CPPL_VERBOSE_REPORT;
  l =0;
  //cap =0;
  array =NULL;
}

//=============================================================================
/*! resize vector */
inline void zcovector::resize(const CPPL_INT& _l)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( _l<0 ){
    ERROR_REPORT;
    std::cerr << "Vector size must be positive integers." << std::endl
              << "Your input was (" << _l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  l =_l;
  delete [] array;
  array =new comple[_l];
}

//=============================================================================
/*! swap two vectors */
inline void swap(zcovector& u, zcovector& v)
{CPPL_VERBOSE_REPORT;
  CPPL_INT u_l =u.l;
  comple* u_array =u.array;
  u.l=v.l; u.array=v.array;
  v.l=u_l; v.array=u_array;
}

//=============================================================================
/*! convert user object to smart-temporary object */
inline _zcovector _(zcovector& vec)
{CPPL_VERBOSE_REPORT;
  _zcovector newvec;
  
  //////// shallow copy ////////
  newvec.l =vec.l;
  newvec.array =vec.array;
  
  //////// nullify ////////
  vec.l =0;
  vec.array =NULL;
  
  return newvec;
}
