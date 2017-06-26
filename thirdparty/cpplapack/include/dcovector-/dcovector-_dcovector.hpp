//=============================================================================
/*! dcovector=_dcovector operator */
inline dcovector& dcovector::operator=(const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
  shallow_copy(vec);
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector+=_dcovector operator */
inline dcovector& dcovector::operator+=(const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l!=vec.l ){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a sumation." << std::endl
              << "Your input was (" << l << ") += (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<l; i++){ array[i]+=vec.array[i]; }
  
  vec.destroy();
  return *this;
}

//=============================================================================
/*! dcovector operator-= */
inline dcovector& dcovector::operator-=(const _dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( l!=vec.l ){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a subtraction." << std::endl
              << "Your input was (" << l << ") -= (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<l; i++){ array[i]-=vec.array[i]; }
  
  vec.destroy();
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector+dcovector operator */
inline _dcovector operator+(const dcovector& vecA, const _dcovector& vecB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vecA.l!=vecB.l){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a sumation." << std::endl
              << "Your input was (" << vecA.l << ") + (" << vecB.l << ")." << std::endl;
    exit(1);
  }
  
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<vecA.l; i++){ vecB.array[i]+=vecA.array[i]; }
  
  return vecB;
}

//=============================================================================
/*! dcovector-dcovector operator */
inline _dcovector operator-(const dcovector& vecA, const _dcovector& vecB)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(vecA.l!=vecB.l){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a subtraction." << std::endl
              << "Your input was (" << vecA.l << ") - (" << vecB.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  for(CPPL_INT i=0; i<vecA.l; i++){
    vecB.array[i] =vecA.array[i]-vecB.array[i];
  }
  
  return vecB;
}

//=============================================================================
/*! dcovector^T*dcovector operator (inner product) */
inline double operator%(const dcovector& vecA, const _dcovector& vecB)
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
  
  double val =ddot_( &vecA.l, vecA.array, &inc, vecB.array, &inc );
  
  vecB.destroy();
  return val;
}
