//=============================================================================
/*! return transposed dgbmatrix */
inline _dgbmatrix t(const dgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dgbmatrix newmat(mat.n, mat.m, mat.ku, mat.kl);
  
  for(CPPL_INT i=0; i<newmat.m; i++){
    const CPPL_INT jmax =std::min(newmat.n,i+newmat.ku+1);
    for(CPPL_INT j=std::max(CPPL_INT(0),i-newmat.kl); j<jmax; j++){
      newmat(i,j) =mat(j,i);
    }
  }
  
  return _(newmat);
}

//=============================================================================
/*! return its inverse matrix */
inline _dgematrix i(const dgbmatrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(mat.m!=mat.n){
    ERROR_REPORT;
    std::cerr << "This matrix is not square and has no inverse matrix." << std::endl
              << "Your input was (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  dgbmatrix mat_cp(mat);
  dgematrix mat_inv(mat.m,mat.n);
  mat_inv.identity();
  mat_cp.dgbsv(mat_inv);
  
  return _(mat_inv);
}
