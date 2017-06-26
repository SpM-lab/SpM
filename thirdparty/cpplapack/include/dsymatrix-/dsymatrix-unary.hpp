//=============================================================================
/*! +dsymatrix operator */
inline const dsymatrix& operator+(const dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -dsymatrix operator */
inline _dsymatrix operator-(const dsymatrix& mat)
{CPPL_VERBOSE_REPORT;
  dsymatrix newmat(mat.n);
  
  const CPPL_INT size =newmat.n*newmat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =-mat.array[i];
  }
  
  return _(newmat);
}
