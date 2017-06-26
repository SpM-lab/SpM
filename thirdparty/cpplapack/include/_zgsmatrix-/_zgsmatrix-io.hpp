//=============================================================================
/*! operator() for const object */
inline comple _zgsmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// search (i,j) component ////////
  const std::vector<CPPL_INT>::iterator rows_i_end =rows[i].end();
  for(std::vector<CPPL_INT>::iterator p=rows[i].begin(); p!=rows_i_end; p++){
    if( data[*p].j==j ){ return data[*p].v; }
  }
  
  //////// (i,j) component was not found ////////
  return comple(0.0,0.0);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const _zgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      std::vector<CPPL_INT>::iterator q;
      const std::vector<CPPL_INT>::iterator mat_rows_i_end =mat.rows[i].end();
      for(q=mat.rows[i].begin(); q!=mat_rows_i_end; q++){
        if( mat.data[*q].j==j ){ break; }
      }
      if(q!=mat_rows_i_end){ s << " " << mat.data[*q].v; }
      else{ s << " x"; }
    }
    s << std::endl;
  }
  
  mat.destroy();
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void _zgsmatrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#zgsmatrix " << m << " " << n << " " << data.size() << std::endl;
  
  const std::vector<zcomponent>::const_iterator data_end =data.end();
  for(std::vector<zcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    ofs << it->i << " " << it->j << " " << it->v << std::endl;
  }
  
  ofs.close();
  destroy();
}
