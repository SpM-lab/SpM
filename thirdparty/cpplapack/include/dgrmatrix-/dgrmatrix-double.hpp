//=============================================================================
/*! dgrmatrix*=double operator */
inline dgrmatrix& dgrmatrix::operator*=(const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<double>::iterator a_end =a.end();
  for(std::vector<double>::iterator it=a.begin(); it!=a_end; it++){
    (*it) *=d;
  }
  
  return *this;
}

//=============================================================================
/*! dgrmatrix/=double operator */
inline dgrmatrix& dgrmatrix::operator/=(const double& d)
{CPPL_VERBOSE_REPORT;
  const std::vector<double>::iterator a_end =a.end();
  for(std::vector<double>::iterator it=a.begin(); it!=a_end; it++){
    (*it) /=d;
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dgrmatrix*double operator */
inline dgrmatrix operator*(const dgrmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgrmatrix newmat =mat;
  
  const std::vector<double>::iterator newmat_a_end =newmat.a.end();
  for(std::vector<double>::iterator it=newmat.a.begin(); it!=newmat_a_end; it++){
    (*it) *=d;
  }
  
  return newmat;
}

//=============================================================================
/*! dgrmatrix/double operator */
inline dgrmatrix operator/(const dgrmatrix& mat, const double& d)
{CPPL_VERBOSE_REPORT;
  dgrmatrix newmat =mat;
  
  const std::vector<double>::iterator newmat_a_end =newmat.a.end();
  for(std::vector<double>::iterator it=newmat.a.begin(); it!=newmat_a_end; it++){
    (*it) /=d;
  }
  
  return newmat;
}
