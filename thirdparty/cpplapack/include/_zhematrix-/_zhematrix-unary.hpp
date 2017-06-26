//=============================================================================
/*! +_zhematrix operator */
inline const _zhematrix& operator+(const _zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -_zhematrix operator */
inline _zhematrix operator-(const _zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT j=0; j<mat.n; j++){
    for(CPPL_INT i=j; i<mat.n; i++){
      mat.darray[j][i] =-mat.darray[j][i];
    }
  }
  
  return mat;
}
