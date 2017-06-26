//=============================================================================
/*! +zgsmatrix operator */
inline const zgsmatrix& operator+(const zgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -zgsmatrix operator */
inline _zgsmatrix operator-(const zgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat(mat);
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v =-it->v;
  }
  
  return _(newmat);
}
