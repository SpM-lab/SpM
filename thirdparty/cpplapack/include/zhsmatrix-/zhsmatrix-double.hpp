//=============================================================================
/*! zhsmatrix*=double operator */
inline zhsmatrix& zhsmatrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator data_end =data.end();
  for(std::vector<zcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v *=d;
  }
  
  return *this;
}

//=============================================================================
/*! zhsmatrix/=double operator */
inline zhsmatrix& zhsmatrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator data_end =data.end();
  for(std::vector<zcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v /=d;
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhsmatrix*double operator */
inline _zhsmatrix operator*(const zhsmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  zhsmatrix newmat(mat);
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v *=d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zhsmatrix/double operator */
inline _zhsmatrix operator/(const zhsmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  zhsmatrix newmat(mat);
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v /=d;
  }
  
  return _(newmat);
}
