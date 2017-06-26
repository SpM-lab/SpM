//=============================================================================
/*! +dgematrix operator */
inline const dgematrix& operator+(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -dgematrix operator */
inline _dgematrix operator-(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(mat.m,mat.n);
  for(CPPL_INT i=0; i<newmat.m*newmat.n; i++){ newmat.array[i]=-mat.array[i]; }
  return _(newmat);
}
