//=============================================================================
/*! return a transposed column vector */
inline _zcovector t(const _zrovector& rovec)
{CPPL_VERBOSE_REPORT;
  _zcovector covec;
  covec.l =rovec.l;
  delete [] covec.array;
  covec.array =rovec.array;
  
  return covec;
}

//=============================================================================
/*! return its conjugated vector */
inline _zrovector conj(const _zrovector& vec)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<vec.l; i++){
    vec(i) =std::conj(vec(i));
  }
  return vec;
}

//=============================================================================
/*! return a conjugate transposed column vector */
inline _zcovector conjt(const _zrovector& rovec)
{CPPL_VERBOSE_REPORT;
  zcovector covec(rovec.l);
  for(CPPL_INT i=0; i<rovec.l; i++){
    covec(i) =std::conj(rovec(i));
  }
  
  rovec.destroy();
  return _(covec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return its Euclidean norm */
inline double nrm2(const _zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  double val =dznrm2_(&vec.l, vec.array, &inc);
  vec.destroy();
  return val;
}

//=============================================================================
/*! return the index of element having the largest absolute value
  in 0-based numbering system */
inline CPPL_INT idamax(const _zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  CPPL_INT i =izamax_(&vec.l, vec.array, &inc) -1;
  vec.destroy();
  return i;
}

//=============================================================================
/*! return its largest absolute value */
inline comple damax(const _zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  comple val =vec.array[izamax_(&vec.l, vec.array, &inc) -1];
  vec.destroy();
  return val;
}
