//=============================================================================
/*! +dssmatrix operator */
inline const dssmatrix& operator+(const dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -dssmatrix operator */
inline _dssmatrix operator-(const dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dssmatrix newmat(mat);
  
  const std::vector<dcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<dcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v =-it->v;
  }
  
  return _(newmat);
}
