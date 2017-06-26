//=============================================================================
/*! cast to _zcovector */
inline _zcovector dcovector::to_zcovector() const
{CPPL_VERBOSE_REPORT;
  zcovector newvec(l);
  
  for(CPPL_INT i=0; i<l; i++){
    newvec.array[i] =comple(array[i], 0.);
  }
  
  return _(newvec);
}
