//=============================================================================
/*! _zhematrix*comple operator */
inline _zgematrix operator*(const _zhematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat =mat.to_zgematrix();
  
  CPPL_INT size =mat.n*mat.n;
  CPPL_INT inc =1;
  zscal_(&size, &d, newmat.array, &inc);
  
  return _(newmat);
}

//=============================================================================
/*! zhematrix/comple operator */
inline _zgematrix operator/(const _zhematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat =mat.to_zgematrix();
  
  CPPL_INT size =mat.n*mat.n;
  comple dinv =1./d;
  CPPL_INT inc =1;
  zscal_(&size, &dinv, newmat.array, &inc);
  
  return _(newmat);
}
