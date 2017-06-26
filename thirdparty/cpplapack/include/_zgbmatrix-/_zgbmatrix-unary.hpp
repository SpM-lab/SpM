//=============================================================================
/*! +_zgbmatrix operator */
inline const _zgbmatrix& operator+(const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -_zgbmatrix operator */
inline _zgbmatrix operator-(const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<(mat.kl+mat.ku+1)*mat.n; i++){
    mat.array[i] =-mat.array[i];
  }
  return mat;
}
