//=============================================================================
/*! return a transposed column vector */
inline _dcovector t(const drovector& rovec)
{CPPL_VERBOSE_REPORT;
  dcovector covec(rovec.l);
  
  CPPL_INT inc =1;
  dcopy_(&rovec.l, rovec.array, &inc, covec.array, &inc);
  
  return _(covec);
}

//=============================================================================
/*! return its Euclidean norm */
inline double nrm2(const drovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return dnrm2_(&vec.l, vec.array, &inc);
}

//=============================================================================
/*! return the index of element having the largest absolute value
  in 0-based numbering system */
inline CPPL_INT idamax(const drovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return idamax_(&vec.l, vec.array, &inc) -1;
}

//=============================================================================
/*! return its largest absolute value */
inline double damax(const drovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return vec.array[idamax_(&vec.l, vec.array, &inc) -1];
}
