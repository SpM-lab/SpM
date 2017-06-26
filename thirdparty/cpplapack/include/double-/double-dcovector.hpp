//=============================================================================
/*! double*dcovector operator */
inline _dcovector operator*(const double& d, const dcovector& vec)
{CPPL_VERBOSE_REPORT;
  dcovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =d*vec.array[i];
  }
  
  return _(newvec);
}
