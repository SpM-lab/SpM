//=============================================================================
/*! double*zrovector operator */
inline _zrovector operator*(const double& d, const zrovector& vec)
{CPPL_VERBOSE_REPORT;
  zrovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =d*vec.array[i];
  }
  
  return _(newvec);
}
