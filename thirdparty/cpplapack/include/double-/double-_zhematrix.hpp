//=============================================================================
/*! double*_zhematrix operator */
inline _zhematrix operator*(const double& d, const _zhematrix& mat)
{CPPL_VERBOSE_REPORT;
  CPPL_INT size =mat.n*mat.n;
  CPPL_INT inc =1;
  zdscal_(&size, &d, mat.array, &inc);
  return mat;
}
