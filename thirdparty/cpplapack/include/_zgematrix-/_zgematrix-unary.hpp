//=============================================================================
/*! +_zgematrix operator */
inline const _zgematrix& operator+(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -_zgematrix operator */
inline _zgematrix operator-(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m*mat.n; i++){ mat.array[i]=-mat.array[i]; }
  return mat;
}
