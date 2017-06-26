//=============================================================================
/*! +_drovector operator */
inline const _drovector& operator+(const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  return vec;
}

//=============================================================================
/*! -_drovector operator */
inline _drovector operator-(const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<vec.l; i++){
    vec.array[i] =-vec.array[i];
  }
  return vec;
}
