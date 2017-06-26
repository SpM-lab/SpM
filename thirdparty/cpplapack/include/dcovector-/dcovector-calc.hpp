//=============================================================================
/*! return a transposed row vector */
inline _drovector t(const dcovector& covec)
{CPPL_VERBOSE_REPORT;
  drovector rovec(covec.l);
  
  CPPL_INT inc =1;
  
  dcopy_(&covec.l, covec.array, &inc, rovec.array, &inc);
  
  return _(rovec);
}

//=============================================================================
/*! return its Euclidean norm */
inline double nrm2(const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return dnrm2_(&vec.l, vec.array, &inc);
}

//=============================================================================
/*! return the index of element having the largest absolute value
 in 0-based numbering system */
inline CPPL_INT idamax(const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc=1;
  return idamax_(&vec.l, vec.array, &inc) -1;
}

//=============================================================================
/*! return its largest absolute value */
inline double damax(const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return vec.array[idamax_(&vec.l, vec.array, &inc) -1];
}
