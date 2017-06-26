//=============================================================================
/*! zgsmatrix*=comple operator */
inline zgsmatrix& zgsmatrix::operator*=(const comple& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<zcomponent>::iterator data_end =data.end();
  for(std::vector<zcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v *=d;
  }
  
  return *this;
}

//=============================================================================
/*! zgsmatrix/=comple operator */
inline zgsmatrix& zgsmatrix::operator/=(const comple& d)
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
/*! zgsmatrix*comple operator */
inline _zgsmatrix operator*(const zgsmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat(mat);
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v *=d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! zgsmatrix/comple operator */
inline _zgsmatrix operator/(const zgsmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat(mat);
  
  const std::vector<zcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<zcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v /=d;
  }
  
  return _(newmat);
}
