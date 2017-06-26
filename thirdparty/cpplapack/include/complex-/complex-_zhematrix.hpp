//=============================================================================
/*! comple*_zhematrix operator */
inline _zgematrix operator*(const comple& d, const _zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  zgematrix newmat =mat.to_zgematrix();
  
  CPPL_INT size =mat.n*mat.n;
  CPPL_INT inc =1;
  zscal_(&size, &d, newmat.array, &inc);
  
  mat.destroy();
  return _(newmat);
}
