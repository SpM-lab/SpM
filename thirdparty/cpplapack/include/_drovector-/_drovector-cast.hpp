//=============================================================================
/*! cast to _zrovector */
inline _zrovector _drovector::to_zrovector() const
{CPPL_VERBOSE_REPORT;
  zrovector newvec(l);
  
  for(CPPL_INT i=0; i<l; i++){
    newvec.array[i] =comple(array[i], 0.);
  }
  
  destroy();
  return _(newvec);
}
