//=============================================================================
/*! +dgsmatrix operator */
inline const dgsmatrix& operator+(const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -dgsmatrix operator */
inline _dgsmatrix operator-(const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dgsmatrix newmat(mat);
  
  const std::vector<dcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<dcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v =-it->v;
  }
  
  return _(newmat);
}
