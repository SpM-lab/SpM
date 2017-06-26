//=============================================================================
/*! _dgbmatrix*double operator */
inline _dgbmatrix operator*(const _dgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  CPPL_INT inc =1;
  dscal_(&size, &d, mat.array, &inc);
  return mat;
}

//=============================================================================
/*! _dgbmatrix/double operator */
inline _dgbmatrix operator/(const _dgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  double dinv =1./d;
  CPPL_INT inc =1;
  dscal_(&size, &dinv, mat.array, &inc);
  return mat;
}
