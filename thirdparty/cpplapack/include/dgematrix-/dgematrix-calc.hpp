//=============================================================================
/*! return transposed dgematrix */
inline _dgematrix t(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(mat.n,mat.m);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =mat(j,i);
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! return its inverse matrix */
inline _dgematrix i(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.m!=mat.n){
    ERROR_REPORT;
    std::cerr << "This matrix is not square and has no inverse matrix." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgematrix mat_cp(mat), mat_inv(mat.m,mat.n);
  mat_inv.identity();
  mat_cp.dgesv(mat_inv);
  
  return _(mat_inv);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! search the index of element having the largest absolute value
  in 0-based numbering system */
inline void idamax(CPPL_INT& i, CPPL_INT& j, const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT mn =mat.m*mat.n;
  CPPL_INT inc =1;
  CPPL_INT index =idamax_(&mn, mat.array, &inc) -1;
  i =index%mat.m;
  j =index/mat.m;
}

//=============================================================================
/*! return its largest absolute value */
inline double damax(const dgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT mn =mat.m*mat.n;
  CPPL_INT inc =1;
  return mat.array[idamax_(&mn, mat.array, &inc) -1];
}
