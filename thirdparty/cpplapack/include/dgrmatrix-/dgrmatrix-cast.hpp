//=============================================================================
/*! cast to _zgrmatrix */
/*
inline _zgrmatrix dgrmatrix::to_zgrmatrix() const
{CPPL_VERBOSE_REPORT;
  zgrmatrix newmat;
  newmat.m =m;
  newmat.n =n;
  newmat.ia =ia;
  newmat.ja =ja;
  
  newmat.a.resize(a.size());
  const size_t a_size =a.size();
  for(size_t k=0; k<a_size; k++){
    newmat.a[k] =comple(a[k],0.);
  }
  
  return _(newmat);
}
*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! convert to _dgematrix */
inline _dgematrix dgrmatrix::to_dgematrix() const
{CPPL_VERBOSE_REPORT;
  dgematrix newmat(m,n);
  newmat.zero();
  
  for(CPPL_INT i=0; i<m; i++){
    int k_beg =ia[i]-1;
    int k_end =ia[i+1]-1;
    for(int k=k_beg; k<k_end; k++){
      int j =ja[k]-1;
      newmat(i,j) =a[k];
    }
  }
  
  return _(newmat);
}
