//=============================================================================
/*! return a transposed column vector */
inline _dcovector t(const _drovector& rovec)
{CPPL_VERBOSE_REPORT;
  _dcovector covec;
  covec.l =rovec.l;
  covec.cap =rovec.cap;
  delete [] covec.array;
  covec.array =rovec.array;
  
  rovec.nullify();
  return covec;
}

//=============================================================================
/*! return its Euclidean norm */
inline double nrm2(const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  double val =dnrm2_(&vec.l, vec.array, &inc);
  vec.destroy();
  return val;
}

//=============================================================================
/*! return the index of element having the largest absolute value
  in 0-based numbering system */
inline CPPL_INT idamax(const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  CPPL_INT i =idamax_(&vec.l, vec.array, &inc) -1;
  vec.destroy();
  return i;
}

//=============================================================================
/*! return its largest absolute value */
inline double damax(const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  double val =vec.array[idamax_(&vec.l, vec.array, &inc) -1];
  vec.destroy();
  return val;
}
