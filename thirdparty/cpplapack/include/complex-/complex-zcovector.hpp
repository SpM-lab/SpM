//=============================================================================
/*! comple*zcovector operator */
inline _zcovector operator*(const comple& d, const zcovector& vec)
{CPPL_VERBOSE_REPORT;
  zcovector newvec(vec.l);
  
  for(CPPL_INT i=0; i<vec.l; i++){
    newvec.array[i] =d*vec.array[i];
  }
  
  return _(newvec);
}
