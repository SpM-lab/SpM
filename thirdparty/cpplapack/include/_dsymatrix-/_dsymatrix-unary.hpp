//=============================================================================
/*! +_dsymatrix operator */
inline const _dsymatrix& operator+(const _dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -_dsymatrix operator */
inline _dsymatrix operator-(const _dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  const CPPL_INT size =mat.n*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    mat.array[i] =-mat.array[i];
  }
  
  return mat;
}
