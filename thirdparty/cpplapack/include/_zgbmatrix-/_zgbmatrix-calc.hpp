//=============================================================================
/*! return its transposed zgbmatrix */
inline _zgbmatrix t(const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(mat.n, mat.m, mat.ku, mat.kl);
  for(CPPL_INT i=0; i<newmat.m; i++){
    const CPPL_INT jmax =std::min(newmat.n,i+newmat.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-newmat.kl); j<jmax; j++){
      newmat(i,j) =mat(j,i);
    }
  }
  
  mat.destroy();
  return _(newmat);
}

//=============================================================================
/*! return its inverse matrix */
inline _zgematrix i(const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.m!=mat.n){
    ERROR_REPORT;
    std::cerr << "This matrix is not square and has no inverse matrix." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  zgbmatrix mat_cp(mat);
  zgematrix mat_inv(mat_cp.m,mat_cp.n);
  mat_inv.identity();
  mat_cp.zgbsv(mat_inv);
  
  return _(mat_inv);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return its conjugate matrix */
inline _zgbmatrix conj(const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m; i++){
    const CPPL_INT jmax =std::min(mat.n,i+mat.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-mat.kl); j<jmax; j++){
      mat(i,j) =std::conj(mat(i,j));
    }
  }
  
  return mat;
}

//=============================================================================
/*! return its conjugate transposed zgbmatrix */
inline _zgbmatrix conjt(const _zgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  zgbmatrix newmat(mat.n, mat.m, mat.ku, mat.kl);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    const CPPL_INT jmax =std::min(newmat.n,i+newmat.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-newmat.kl); j<jmax; j++){
      newmat(i,j) =std::conj(mat(j,i));
    }
  }
  
  mat.destroy();
  return _(newmat);
}
