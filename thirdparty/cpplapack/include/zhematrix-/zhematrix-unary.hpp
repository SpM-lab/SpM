//=============================================================================
/*! +zhematrix operator */
inline const zhematrix& operator+(const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -zgematrix operator */
inline _zhematrix operator-(const zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(mat.n);
  
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      newmat(i,j) =-mat(i,j);
    }
  }
  
  return _(newmat);
}
