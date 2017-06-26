//=============================================================================
/*! zrovector=zrovector operator */
inline zrovector& zrovector::operator=(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  if(array!=vec.array){ // if it is NOT self substitution
    copy(vec);
  }
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector+=zrovector operator */
inline zrovector& zrovector::operator+=(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l!=vec.l ){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a sumation." << std::endl
              << "Your input was (" << l << ") += (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<l; i++){
    array[i] += vec.array[i];
  }
  
  return *this;
}

//=============================================================================
/*! zrovector operator-= */
inline zrovector& zrovector::operator-=(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l!=vec.l ){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a subtraction." << std::endl
              << "Your input was (" << l << ") -= (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<l; i++){
    array[i] -= vec.array[i];
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector+zrovector operator */
inline _zrovector operator+(const zrovector& vecA, const zrovector& vecB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vecA.l!=vecB.l){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a sumation." << std::endl
              << "Your input was (" << vecA.l << ") + (" << vecB.l << ")." << std::endl;
    exit(1);
  }
  
#endif//CPPL_DEBUG
  
  zrovector newvec(vecA.l);
  
  for(CPPL_INT i=0; i<newvec.l; i++){
    newvec.array[i] =vecA.array[i]+vecB.array[i];
  }
  
  return _(newvec);
}

//=============================================================================
/*! zrovector-zrovector operator */
inline _zrovector operator-(const zrovector& vecA, const zrovector& vecB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vecA.l!=vecB.l){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a subtraction." << std::endl
              << "Your input was (" << vecA.l << ") - (" << vecB.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zrovector newvec(vecA.l);
  
  for(CPPL_INT i=0; i<newvec.l; i++){
    newvec.array[i] =vecA.array[i]-vecB.array[i];
  }
  
  return _(newvec);
}

//=============================================================================
/*! zrovector^T*zrovector operator (inner product) */
inline comple operator%(const zrovector& vecA, const zrovector& vecB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vecA.l!=vecB.l){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a dot product." << std::endl
              << "Your input was (" << vecA.l << ") % (" << vecB.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  CPPL_INT inc =1;
  comple val =zdotu_( &vecA.l, vecA.array, &inc, vecB.array, &inc );
  return val;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return Hadamerd product */
inline _zrovector hadamerd(const zrovector& vecA, const zrovector& vecB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( vecA.l!=vecB.l ){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make Hadamerd product." << std::endl
              << "Your input was (" << vecA.l << ") and (" << vecB.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  zrovector newvec(vecA.l);
  for(CPPL_INT i=0; i<newvec.l; i++){
    newvec(i) =vecA(i)*vecB(i);
  }
  return _(newvec);
}
