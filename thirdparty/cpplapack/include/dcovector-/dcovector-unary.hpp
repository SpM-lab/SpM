//=============================================================================
/*! +dcovector operator */
inline const dcovector& operator+(const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  return vec;
}

//=============================================================================
/*! -dcovector operator */
inline _dcovector operator-(const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  dcovector newvec(vec.l);
  for(CPPL_INT i=0; i<newvec.l; i++){ newvec.array[i]=-vec.array[i]; }
  
  return _(newvec);
}
