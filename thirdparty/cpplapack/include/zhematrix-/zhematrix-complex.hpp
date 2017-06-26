//=============================================================================
/*! zhematrix*comple operator */
inline _zgematrix operator*(const zhematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  mat.complete();
  zgematrix newmat(mat.n, mat.n);
  
  const CPPL_INT size =mat.n*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]*d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zhematrix/comple operator */
inline _zgematrix operator/(const zhematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  mat.complete();
  zgematrix newmat(mat.n, mat.n);

  const CPPL_INT size =mat.n*mat.n;
  for(CPPL_INT i=0; i<size; i++){
    newmat.array[i] =mat.array[i]/d;
  }
  
  return _(newmat);
}
