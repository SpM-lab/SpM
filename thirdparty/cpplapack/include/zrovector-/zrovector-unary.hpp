//=============================================================================
/*! +zrovector operator */
inline const zrovector& operator+(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  return vec;
}

//=============================================================================
/*! -zrovector operator */
inline _zrovector operator-(const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  zrovector newvec(vec.l);
  for(CPPL_INT i=0; i<newvec.l; i++){ newvec.array[i]=-vec.array[i]; }
  
  return _(newvec);
}
