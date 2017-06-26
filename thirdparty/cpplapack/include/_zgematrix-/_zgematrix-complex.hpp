//=============================================================================
/*! _zgematrix*comple operator */
inline _zgematrix operator*(const _zgematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  CPPL_INT inc =1;
  zscal_(&size, &d, mat.array, &inc);
  return mat;
}

//=============================================================================
/*! _zgematrix/comple operator */
inline _zgematrix operator/(const _zgematrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.m*mat.n;
  comple dinv =1./d;
  CPPL_INT inc =1;
  zscal_(&size, &dinv, mat.array, &inc);
  return mat;
}
