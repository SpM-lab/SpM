//=============================================================================
/*! _zgbmatrix*double operator */
inline _zgbmatrix operator*(const _zgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, mat.array, &inc);
  return mat;
}

//=============================================================================
/*! _zgbmatrix/double operator */
inline _zgbmatrix operator/(const _zgbmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =(mat.kl+mat.ku+1)*mat.n;
  double dinv =1./d;
  CPPL_INT inc =1;
  zdscal_(&size, &dinv, mat.array, &inc);
  return mat;
}
