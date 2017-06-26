//=============================================================================
/*! +_zcovector operator */
inline const _zcovector& operator+(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  return vec;
}

//=============================================================================
/*! -_zcovector operator */
inline _zcovector operator-(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<vec.l; i++){ vec.array[i]=-vec.array[i]; }
  return vec;
}
