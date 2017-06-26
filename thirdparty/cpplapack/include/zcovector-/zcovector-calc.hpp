//=============================================================================
/*! return a transposed row vector */
inline _zrovector t(const zcovector& covec)
{CPPL_VERBOSE_REPORT;
  zrovector rovec(covec.l);
  
  CPPL_INT inc =1;
  zcopy_(&covec.l, covec.array, &inc, rovec.array, &inc);
  
  return _(rovec);
}
//=============================================================================
/*! return its conjugated vector */
inline _zcovector conj(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  zcovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec(i) =std::conj(vec(i));
  }
  
  return _(newvec);
}

//=============================================================================
/*! return a conjugate transposed row vector */
inline _zrovector conjt(const zcovector& covec)
{CPPL_VERBOSE_REPORT;
  zrovector rovec(covec.l);
  
  for(CPPL_INT i=0; i<covec.l; i++){
    rovec(i) =std::conj(covec(i));
  }
  
  return _(rovec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return its Euclidean norm */
inline double nrm2(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return dznrm2_(&vec.l, vec.array, &inc);
}

//=============================================================================
/*! return the index of element having the largest absolute value
 in 0-based numbering system */
inline CPPL_INT idamax(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return izamax_(&vec.l, vec.array, &inc) -1;
}

//=============================================================================
/*! return its largest absolute value */
inline comple damax(const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  CPPL_INT inc =1;
  return vec.array[izamax_(&vec.l, vec.array, &inc) -1];
}
