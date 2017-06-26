//=============================================================================
/*! return transposed zgematrix */
inline _zgematrix t(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.n,mat.m);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =mat(j,i);
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! return its inverse matrix */
inline _zgematrix i(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.m!=mat.n){
    ERROR_REPORT;
    std::cerr << "This matrix is not square and has no inverse matrix." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix mat_cp(mat), mat_inv(mat.m,mat.n);
  mat_inv.identity();
  mat_cp.zgesv(mat_inv);
  
  return _(mat_inv);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return its conjugate matrix */
inline _zgematrix conj(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.m,mat.n);
  
  for(CPPL_INT i=0; i<mat.m; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      newmat(i,j) =std::conj(mat(i,j));
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! return its conjugate transposed matrix */
inline _zgematrix conjt(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.n,mat.m);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =std::conj(mat(j,i));
    }
  }
  
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! search the index of element having the largest absolute value
  in 0-based numbering system */
inline void idamax(CPPL_INT& i, CPPL_INT& j, const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  CPPL_INT index =izamax_(&size, mat.array, &inc) -1;
  i =index%mat.m;
  j =index/mat.m;
}

//=============================================================================
/*! return its largest absolute value */
inline comple damax(const zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  return mat.array[izamax_(&size, mat.array, &inc) -1];
}
