//=============================================================================
/*! +dgrmatrix operator */
inline const dgrmatrix& operator+(const dgrmatrix& mat)
{CPPL_VERBOSE_REPORT;
  return mat;
}

//=============================================================================
/*! -dgrmatrix operator */
inline dgrmatrix operator-(const dgrmatrix& mat)
{CPPL_VERBOSE_REPORT;
  dgrmatrix newmat =mat;
  
  const std::vector<double>::iterator newmat_a_end =newmat.a.end();
  for(std::vector<double>::iterator it=newmat.a.begin(); it!=newmat_a_end; it++){
    (*it) =-(*it);
  }
  
  return newmat;
}
