//=============================================================================
/*! _zgbmatrix*comple operator */
inline _zgbmatrix operator*(const _zgbmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  CPPL_INT inc =1;
  zscal_(&size, &d, mat.array, &inc);
  return mat;
}

//=============================================================================
/*! _zgbmatrix/comple operator */
inline _zgbmatrix operator/(const _zgbmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  comple dinv =1./d;
  CPPL_INT inc =1;
  zscal_(&size, &dinv, mat.array, &inc);
  return mat;
}
