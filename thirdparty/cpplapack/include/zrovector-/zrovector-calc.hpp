//=============================================================================
/*! return a transposed column vector */
inline _zcovector t(const zrovector& rovec)
{CPPL_VERBOSE_REPORT;
  zcovector covec(rovec.l);
  CPPL_INT inc =1;
  
  zcopy_(&rovec.l, rovec.array, &inc, covec.array, &inc);
  
  return _(covec);
}

//=============================================================================
/*! return its conjugated vector */
inline _zrovector conj(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  zrovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec(i) =std::conj(vec(i));
  }
  
  return _(newvec);
}

//=============================================================================
/*! return a conjugate transposed column vector */
inline _zcovector conjt(const zrovector& rovec)
{CPPL_VERBOSE_REPORT;
  zcovector covec(rovec.l);
  
  for(CPPL_INT i=0; i<rovec.l; i++){
    covec(i) =std::conj(rovec(i));
  }
  
  return _(covec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return its Euclidean norm */
inline double nrm2(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return dznrm2_(&vec.l, vec.array, &inc);
}

//=============================================================================
/*! return the index of element having the largest absolute value
  in 0-based numbering system */
inline CPPL_INT idamax(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return izamax_(&vec.l, vec.array, &inc) -1;
}

//=============================================================================
/*! return its largest absolute value */
inline comple damax(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return vec.array[izamax_(&vec.l, vec.array, &inc) -1];
}
