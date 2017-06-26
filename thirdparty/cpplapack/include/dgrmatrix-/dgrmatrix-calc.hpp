//=============================================================================
/*! return its largest absolute value */
inline double damax(const dgrmatrix& mat)
{CPPL_VERBOSE_REPORT;
  //////// exception ////////
  if(mat.a.size()==0){
    return 0.;
  }
  
  //////// find ////////
  const size_t mat_a_size =mat.a.size();
  double amax =0.;
  double vmax;
  for(size_t k=0; k<mat_a_size; k++){
    if( amax < fabs(mat.a[k]) ){
      amax =fabs(mat.a[k]);
      vmax =mat.a[k];
    }
  }
  
  return vmax;
}
