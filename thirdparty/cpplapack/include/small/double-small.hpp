//=============================================================================
/*! double*dcovector_small operator */
template<CPPL_INT l>
inline dcovector_small<l> operator*(const double& v, const dcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  dcovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =v*A(i);
  }
  return X;
}

//=============================================================================
/*! double*drovector_small operator */
template<CPPL_INT l>
inline drovector_small<l> operator*(const double& v, const drovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =v*A(i);
  }
  return X;
}

//=============================================================================
/*! double*dgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline dgematrix_small<m,n> operator*(const double& v, const dgematrix_small<m,n>& A)
{CPPL_VERBOSE_REPORT;
  dgematrix_small<m,n> C;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i,j) =v*A(i,j);
    }
  }
  return C;
}

//=============================================================================
/*! double*dsymatrix_small operator */
template<CPPL_INT n>
inline dsymatrix_small<n> operator*(const double& v, const dsymatrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  dsymatrix_small<n> X;
  for(CPPL_INT k=0; k<(n*(n+1))/2; k++){
    X.array[k] =v*A.array[k];
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! double*zcovector_small operator */
template<CPPL_INT l>
inline zcovector_small<l> operator*(const double& v, const zcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  zcovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =v*A(i);
  }
  return X;
}

//=============================================================================
/*! double*zrovector_small operator */
template<CPPL_INT l>
inline zrovector_small<l> operator*(const double& v, const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =v*A(i);
  }
  return X;
}

//=============================================================================
/*! double*zgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator*(const double& v, const zgematrix_small<m,n>& A)
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
/*! double*zhematrix_small operator */
template<CPPL_INT n>
inline zhematrix_small<n> operator*(const double& v, const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> X;
  for(CPPL_INT k=0; k<(n*(n+1))/2; k++){
    X.array[k] =v*A.array[k];
  }
  return X;
}
