//=============================================================================
/*! return a transposed row vector */
inline _zrovector t(const _zcovector& covec)
{CPPL_VERBOSE_REPORT;
  _zrovector rovec;
  rovec.l =covec.l;
  delete [] rovec.array;
  rovec.array =covec.array;
  
  return rovec;
}

//=============================================================================
/*! return its conjugated vector */
inline _zcovector conj(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<vec.l; i++){
    vec(i) =std::conj(vec(i));
  }
  return vec;
}

//=============================================================================
/*! return a conjugate transposed row vector */
inline _zrovector conjt(const _zcovector& covec)
{CPPL_VERBOSE_REPORT;
  zrovector rovec(covec.l);
  
  for(CPPL_INT i=0; i<covec.l; i++){
    rovec(i) =std::conj(covec(i));
  }
  
  covec.destroy();
  return _(rovec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return its Euclidean norm */
inline double nrm2(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  double val =dznrm2_(&vec.l, vec.array, &inc);
  vec.destroy();
  return val;
}

//=============================================================================
/*! return the index of element having the largest absolute value
 in 0-based numbering system */
inline CPPL_INT idamax(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  CPPL_INT i =izamax_(&vec.l, vec.array, &inc) -1;
  vec.destroy();
  return i;
}

//=============================================================================
/*! return its largest absolute value */
inline comple damax(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  comple val =vec.array[izamax_(&vec.l, vec.array, &inc) -1];
  vec.destroy();
  return val;
}
