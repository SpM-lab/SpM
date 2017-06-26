//=============================================================================
/*! +zgematrix operator */
inline const zgematrix& operator+(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -zgematrix operator */
inline _zgematrix operator-(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.m,mat.n);
  for(CPPL_INT i=0; i<newmat.m*newmat.n; i++){ newmat.array[i]=-mat.array[i]; }
  
  return _(newmat);
}
