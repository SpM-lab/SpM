//=============================================================================
/*! return transposed zgematrix */
inline _zgematrix t(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.n,mat.m);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =mat(j,i);
    }
  }
  
  mat.destroy();
  return _(newmat);
}

//=============================================================================
/*! return its inverse matrix */
inline _zgematrix i(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.m!=mat.n){
    ERROR_REPORT;
    std::cerr << "This matrix is not square and has no inverse matrix." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgematrix mat_cp(mat);
  zgematrix mat_inv(mat_cp.m,mat_cp.n);
  mat_inv.identity();
  mat_cp.zgesv(mat_inv);
  
  return _(mat_inv);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return its conjugate matrix */
inline _zgematrix conj(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      mat(i,j) =std::conj(mat(i,j));
    }
  }
  
  return mat;
}

//=============================================================================
/*! return its conjugate transposed matrix */
inline _zgematrix conjt(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat(mat.n,mat.m);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    for(CPPL_INT j=0; j<newmat.n; j++){
      newmat(i,j) =std::conj(mat(j,i));
    }
  }
  
  mat.destroy();
  return _(newmat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! search the index of element having the largest absolute value
  in 0-based numbering system */
inline void idamax(CPPL_INT& i, CPPL_INT& j, const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  CPPL_INT index =izamax_(&size, mat.array, &inc) -1;
  i =index%mat.m;
  j =index/mat.m;
  
  mat.destroy();
}

//=============================================================================
/*! return its largest absolute value */
inline comple damax(const _zgematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  comple val =mat.array[izamax_(&size, mat.array, &inc) -1];
  
  mat.destroy();
  return val;
}
