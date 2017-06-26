//=============================================================================
/*! double*drovector operator */
inline _drovector operator*(const double& d, const drovector& vec)
{CPPL_VERBOSE_REPORT;
  drovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =d*vec.array[i];
  }
  
  return _(newvec);
}
