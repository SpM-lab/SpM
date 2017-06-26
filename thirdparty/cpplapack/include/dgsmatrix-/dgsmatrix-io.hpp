//=============================================================================
/*! operator() for const object */
inline double dgsmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// search (i,j) component ////
  const std::vector<CPPL_INT>::const_iterator rows_i_end =rows[i].end();
  for(std::vector<CPPL_INT>::const_iterator p=rows[i].begin(); p!=rows_i_end; p++){
    if(data[*p].j==j){ return data[*p].v; }
  }
  
  //// (i,j) component was not found ////
  return 0.0;
}

//=============================================================================
/*! operator() for const object */
inline double& dgsmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j)
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
    if(data[*p].j==j){ return data[*p].v; }
  }
  
  //////// (i,j) component not found ////////
  rows[i].push_back(CPPL_INT(data.size()));
  cols[j].push_back(CPPL_INT(data.size()));
  data.push_back(dcomponent(i,j,0.));
  return data.back().v;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! put value with volume cheack without isListed check */
inline dgsmatrix& dgsmatrix::put(const CPPL_INT& i, const CPPL_INT& j, const double& v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    if( it->i==i && it->j==j ){
      ERROR_REPORT;
      std::cerr << "The required component is already listed." << std::endl
                << "Your input was (" << i << "," << j << "," << v << ")." << std::endl;
      exit(1);
    }
  }
#endif//CPPL_DEBUG
  
  //// push ////
  rows[i].push_back(CPPL_INT(data.size()));
  cols[j].push_back(CPPL_INT(data.size()));
  data.push_back(dcomponent(i,j,v));
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! delete the entry of a component */
inline dgsmatrix& dgsmatrix::del(const CPPL_INT i, const CPPL_INT j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// search (i,j) component ////
  const std::vector<CPPL_INT>::iterator rows_i_end =rows[i].end();
  for(std::vector<CPPL_INT>::iterator p=rows[i].begin(); p!=rows_i_end; p++){
    if(data[*p].j==j){//exists
      //// save position ////
      CPPL_INT c =*p;
      CPPL_INT C =CPPL_INT(data.size()-1);
      
      //// data translation ////
      data[c]=data.back();
      data.pop_back();
      
      //// remove from List ////
      rows[i].erase(p);
      const std::vector<CPPL_INT>::iterator cols_j_end =cols[j].end();
      for(std::vector<CPPL_INT>::iterator q=cols[j].begin(); q!=cols_j_end; q++){
        if(*q==c){ cols[j].erase(q); break; }
      }
      
      //// modify the entry of translated component ////
      CPPL_INT I(data[c].i), J(data[c].j);
      const std::vector<CPPL_INT>::iterator rows_I_end =rows[I].end();
      for(std::vector<CPPL_INT>::iterator q=rows[I].begin(); q!=rows_I_end; q++){
        if(*q==C){ *q=c; break; }
      }
      const std::vector<CPPL_INT>::iterator cols_J_end =cols[J].end();
      for(std::vector<CPPL_INT>::iterator q=cols[J].begin(); q!=cols_J_end; q++){
        if(*q==C){ *q=c; break; }
      }
      return *this;
    }
  }
  
#ifdef  CPPL_DEBUG
  std::cerr << "# [NOTE]@dgsmatrix::del(CPPL_INT&, CPPL_INT&): The required component was not listed. Your input was (" << i << "," << j << ")." << std::endl;
#endif//CPPL_DEBUG
  
  return *this;
}

//=============================================================================
/*! delete the entry of an element */
inline dgsmatrix& dgsmatrix::del(const CPPL_INT c)
{CPPL_VERBOSE_REPORT;
  if( c==CPPL_INT(data.size()-1) ){//if c is the last element
    CPPL_INT i(data[c].i), j(data[c].j);
    const std::vector<CPPL_INT>::iterator rows_i_end =rows[i].end();
    for(std::vector<CPPL_INT>::iterator q=rows[i].begin(); q!=rows_i_end; q++){
      if(*q==c){ rows[i].erase(q); break; }
    }
    const std::vector<CPPL_INT>::iterator cols_j_end =cols[j].end();
    for(std::vector<CPPL_INT>::iterator q=cols[j].begin(); q!=cols_j_end; q++){
      if(*q==c){ cols[j].erase(q); break; }
    }
    data.pop_back();
  }
  
  else{//if c is NOT the last element
    //// data translation ////
    CPPL_INT C =CPPL_INT(data.size()-1);
    CPPL_INT i(data[c].i), j(data[c].j), I(data.back().i), J(data.back().j);
    data[c]=data.back();
    //// remove entry of component ////
    const std::vector<CPPL_INT>::iterator rows_i_end =rows[i].end();
    for(std::vector<CPPL_INT>::iterator q=rows[i].begin(); q!=rows_i_end; q++){
      if(*q==c){ rows[i].erase(q); break; }
    }
    const std::vector<CPPL_INT>::iterator cols_j_end =cols[j].end();
    for(std::vector<CPPL_INT>::iterator q=cols[j].begin(); q!=cols_j_end; q++){
      if(*q==c){ cols[j].erase(q); break; }
    }
    //// modify the entry of translated component ////
    const std::vector<CPPL_INT>::iterator rows_I_end =rows[I].end();
    for(std::vector<CPPL_INT>::iterator q=rows[I].begin(); q!=rows_I_end; q++){
      if(*q==C){ *q=c; break; }
    }
    const std::vector<CPPL_INT>::iterator cols_J_end =cols[J].end();
    for(std::vector<CPPL_INT>::iterator q=cols[J].begin(); q!=cols_J_end; q++){
      if(*q==C){ *q=c; break; }
    }
    //// pop_back ////
    data.pop_back();
  }
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const dgsmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.m; i++){
    const std::vector<CPPL_INT>::const_iterator mat_rows_i_end =mat.rows[i].end();
    for(CPPL_INT j=0; j<mat.n; j++){
      std::vector<CPPL_INT>::const_iterator q;
      for(q=mat.rows[i].begin(); q!=mat_rows_i_end; q++){
        if(mat.data[*q].j==j){ break; }
      }
      if(q!=mat_rows_i_end){ s << " " << mat.data[*q].v; }
      else{ s << " x"; }
    }
    s << std::endl;
  }
  
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void dgsmatrix::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#dgsmatrix " << m << " " << n << " " << data.size() << std::endl;
  
  const std::vector<dcomponent>::const_iterator data_end =data.end();
  for(std::vector<dcomponent>::const_iterator it=data.begin(); it!=data_end; it++){
    ofs << it->i << " " << it->j << " " << it->v << std::endl;
  }
  
  ofs.close();
}

//=============================================================================
inline void dgsmatrix::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }
  
  //////// id ////////
  std::string id;
  s >> id;
  if( id != "dgsmatrix" && id != "#dgsmatrix" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not dgsmatrix." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  //////// m, n, size ////////
  size_t size;
  s >> m >> n >> size;
  resize(m, n);
  data.resize(size);
  
  //////// i, j, v ////////
  CPPL_INT i, j;
  double v;
  for(size_t k=0; k<size; k++){
    s >> i >> j >> v;
    put(i,j, v);
  }
  
  //////// check ////////
  s >> i;
  if(!s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << " ." << std::endl
              << "Most likely, there are too many data components over the context." << std::endl;
    exit(1);
  }
  s.close();
}
