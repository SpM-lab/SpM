//=============================================================================
/*! operator complex*zcovector_small */
template<CPPL_INT l>
inline zcovector_small<l> operator*(const comple& v, const zcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  zcovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =v*A(i);
  }
  return X;
}

//=============================================================================
/*! operator complex*zrovector_small */
template<CPPL_INT l>
inline zrovector_small<l> operator*(const comple& v, const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =v*A(i);
  }
  return X;
}

//=============================================================================
/*! operator complex*zgematrix_small */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator*(const comple& v, const zgematrix_small<m,n>& A)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i,j) =v*A(i,j);
    }
  }
  return C;
}

//=============================================================================
/*! operator complex*zhematrix_small */
template<CPPL_INT n>
inline zhematrix_small<n> operator*(const comple& v, const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> X;
  for(CPPL_INT k=0; k<(n*(n+1))/2; k++){
    X.array[k] =v*A.array[k];
  }
  return X;
}
