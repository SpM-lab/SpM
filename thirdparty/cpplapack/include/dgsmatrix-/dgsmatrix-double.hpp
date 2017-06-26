//=============================================================================
/*! dgsmatrix*=double operator */
inline dgsmatrix& dgsmatrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator data_end =data.end();
  for(std::vector<dcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v *=d;
  }
  
  return *this;
}

//=============================================================================
/*! dgsmatrix/=double operator */
inline dgsmatrix& dgsmatrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<dcomponent>::iterator data_end =data.end();
  for(std::vector<dcomponent>::iterator it=data.begin(); it!=data_end; it++){
    it->v /=d;
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgsmatrix*double operator */
inline _dgsmatrix operator*(const dgsmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgsmatrix newmat(mat);
  
  const std::vector<dcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<dcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v *=d;
  }
  
  return _(newmat);
}

//=============================================================================
/*! dgsmatrix/double operator */
inline _dgsmatrix operator/(const dgsmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgsmatrix newmat(mat);
  
  const std::vector<dcomponent>::iterator newmat_data_end =newmat.data.end();
  for(std::vector<dcomponent>::iterator it=newmat.data.begin(); it!=newmat_data_end; it++){
    it->v /=d;
  }
  
  return _(newmat);
}
