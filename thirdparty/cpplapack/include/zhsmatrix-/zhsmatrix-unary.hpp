//=============================================================================
/*! +zhsmatrix operator */
inline const zhsmatrix& operator+(const zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -zhsmatrix operator */
inline _zhsmatrix operator-(const zhsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  zhsmatrix newmat(mat);
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v =-it->v;
  }
  
  return _(newmat);
}
