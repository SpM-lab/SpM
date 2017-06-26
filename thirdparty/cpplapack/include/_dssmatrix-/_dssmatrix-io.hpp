//=============================================================================
/*! operator() for const object */
inline double _dssmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is (" << n << "," << n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// search (i,j) component ////
  const CPPL_INT ii(std::max(i,j)), jj(std::min(i,j));
  
  const std::vector<CPPL_INT>::iterator line_ii_end =line[ii].end();
  for(std::vector<CPPL_INT>::iterator p=line[ii].begin(); p!=line_ii_end; p++){
    if(data[*p].i==ii && data[*p].j==jj){ return data[*p].v; }
  }
  
  //// (i,j) component was not found ////
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const _dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.n; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      if( i >= j ){
        std::vector<CPPL_INT>::iterator q;
        const std::vector<CPPL_INT>::iterator mat_line_i_end =mat.line[i].end();
        for(q=mat.line[i].begin(); q!=mat_line_i_end; q++){
          if(mat.data[*q].j==j){ break; }
        }
        if(q!=mat_line_i_end){ s << " " << mat.data[*q].v << " "; }
        else{ s << " x "; }
      }
      else{//i<j
        std::vector<CPPL_INT>::iterator q;
        const std::vector<CPPL_INT>::iterator mat_line_i_end =mat.line[i].end();
        for(q=mat.line[i].begin(); q!=mat_line_i_end; q++){
          if(mat.data[*q].j==j){ break; }
        }
        if(q!=mat_line_i_end){ s << "{" << mat.data[*q].v << "}"; }
        else{ s << "{x}"; }
      }
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
inline void _dssmatrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#dssmatrix " << n << " " << data.size() << std::endl;
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    ofs << it->i << " " << it->j << " " << it->v << std::endl;
  }
  
  ofs.close();
  destroy();
}
