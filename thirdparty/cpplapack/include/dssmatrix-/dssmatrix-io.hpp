//=============================================================================
/*! operator() for const object */
inline double dssmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// search (i,j) component ////
  const CPPL_INT ii(std::max(i,j)), jj(std::min(i,j));
  
  const std::vector<CPPL_INT>::const_iterator line_ii_end =line[ii].end();
  for(std::vector<CPPL_INT>::const_iterator p=line[ii].begin(); p!=line_ii_end; p++){
    if(data[*p].i==ii && data[*p].j==jj){
      return data[*p].v;
    }
  }
  
  //// (i,j) component was not found ////
  return 0.0;
}

//=============================================================================
/*! operator() for const object */
inline double& dssmatrix::operator()(const CPPL_INT& i, const CPPL_INT& j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //////// search (i,j) component ////////
  const CPPL_INT ii(std::max(i,j)), jj(std::min(i,j));
  
  const std::vector<CPPL_INT>::iterator line_ii_end =line[ii].end();
  for(std::vector<CPPL_INT>::iterator p=line[ii].begin(); p!=line_ii_end; p++){
    if(data[*p].i==ii && data[*p].j==jj){
      return data[*p].v;
    }
  }
  
  //////// (i,j) component not found ////////
  line[ii].push_back(CPPL_INT(data.size()));
  if(i!=j){//off-diagonal
    line[jj].push_back(CPPL_INT(data.size()));
  }
  data.push_back(dcomponent(ii,jj,0.));
  return data.back().v;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! put value with volume cheack without isListed check */
inline dssmatrix& dssmatrix::put(const CPPL_INT& i, const CPPL_INT& j, const double& v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
  
  if( isListed(i,j) ){
    ERROR_REPORT;
    std::cerr << "The required component is already listed." << std::endl
              << "Your input was (" << i << "," << j << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //// push ////
  const CPPL_INT ii(std::max(i,j)), jj(std::min(i,j));
  line[ii].push_back(CPPL_INT(data.size()));
  if(i!=j){//off-diagonal
    line[jj].push_back(CPPL_INT(data.size()));
  }
  data.push_back(dcomponent(ii,jj,v));
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! delete the entry of a component */
inline dssmatrix& dssmatrix::del(const CPPL_INT i, const CPPL_INT j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || n<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << n << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG

  const CPPL_INT ii(std::max(i,j)), jj(std::min(i,j));
  
  //////// search (i,j) component ////////
  const std::vector<CPPL_INT>::iterator line_ii_end =line[ii].end();
  for(std::vector<CPPL_INT>::iterator p=line[ii].begin(); p!=line_ii_end; p++){
    if(data[*p].i==ii && data[*p].j==jj){//exists
      //// save position ////
      CPPL_INT c =*p;
      CPPL_INT C =CPPL_INT(data.size()-1);
      
      //// data translation ////
      data[c]=data.back();
      data.pop_back();
      
      //// remove from List ////
      line[ii].erase(p);
      if(i!=j){//off-diagonal
        const std::vector<CPPL_INT>::iterator line_jj_end =line[jj].end();
        for(std::vector<CPPL_INT>::iterator pj=line[jj].begin(); pj!=line_jj_end; pj++){
          if(*pj==c){ line[jj].erase(pj); break; }
        }
      }
      
      //// modify the entry of translated component ////
      CPPL_INT I(data[c].i), J(data[c].j);
      const std::vector<CPPL_INT>::iterator line_I_end =line[I].end();
      for(std::vector<CPPL_INT>::iterator q=line[I].begin(); q!=line_I_end; q++){
        if(*q==C){ *q=c; break; }
      }
      if(I!=J){//off-diagonal
        const std::vector<CPPL_INT>::iterator line_J_end =line[J].end();
        for(std::vector<CPPL_INT>::iterator q=line[J].begin(); q!=line_J_end; q++){
          if(*q==C){ *q=c; break; }
        }
      }
      return *this;
    }
  }
  
#ifdef  CPPL_DEBUG
  std::cerr << "# [NOTE]@dssmatrix::del(CPPL_INT&, CPPL_INT&): The required component was not listed. Your input was (" << i << "," << j << ")." << std::endl;
#endif//CPPL_DEBUG
  
  return *this;
}

//=============================================================================
/*! delete the entry of an element */
inline dssmatrix& dssmatrix::del(const CPPL_INT c)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( c<0 || c>=CPPL_INT(data.size()) ){
    ERROR_REPORT;
    std::cerr << "The required element is out of the matrix volume." << std::endl
              << "Your input was (" << c << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  if( c==CPPL_INT(data.size()-1) ){//if c is the last element
    CPPL_INT i(data[c].i), j(data[c].j);
    const std::vector<CPPL_INT>::iterator line_i_end =line[i].end();
    for(std::vector<CPPL_INT>::iterator q=line[i].begin(); q!=line_i_end; q++){
      if(*q==c){ line[i].erase(q); break; }
    }
    if(i!=j){//off-diagonal
      const std::vector<CPPL_INT>::iterator line_j_end =line[j].end();
      for(std::vector<CPPL_INT>::iterator q=line[j].begin(); q!=line_j_end; q++){
        if(*q==c){ line[j].erase(q); break; }
      }
    }
    data.pop_back();
  }
  
  else{//c is NOT the last element
    //// data translation ////
    CPPL_INT C =CPPL_INT(data.size()-1);
    CPPL_INT i(data[c].i), j(data[c].j), I(data.back().i), J(data.back().j);
    data[c]=data.back();
    //std::cerr << "c=" << c   << " i=" << i << " j=" << j << " C=" << vol << " I=" << I << " J=" << J << std::endl;
    
    //// remove entry of component ////
    const std::vector<CPPL_INT>::iterator line_i_end =line[i].end();
    for(std::vector<CPPL_INT>::iterator q=line[i].begin(); q!=line_i_end; q++){
      if(*q==c){ line[i].erase(q); break; }
    }
    if(i!=j){//off-diagonal
      const std::vector<CPPL_INT>::iterator line_j_end =line[j].end();
      for(std::vector<CPPL_INT>::iterator q=line[j].begin(); q!=line_j_end; q++){
        if(*q==c){ line[j].erase(q); break; }
      }
    }
    
    //// modify the entry of translated component ////
    const std::vector<CPPL_INT>::iterator line_I_end =line[I].end();
    for(std::vector<CPPL_INT>::iterator q=line[I].begin(); q!=line_I_end; q++){
      if(*q==C){ *q=c; break; }
    }
    if(I!=J){//off-diagonal
      const std::vector<CPPL_INT>::iterator line_J_end =line[J].end();
      for(std::vector<CPPL_INT>::iterator q=line[J].begin(); q!=line_J_end; q++){
        if(*q==C){ *q=c; break; }
      }
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
inline std::ostream& operator<<(std::ostream& s, const dssmatrix& mat)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<mat.n; i++){
    for(CPPL_INT j=0; j<mat.n; j++){
      if( i >= j ){
        CPPL_INT c =mat.number(i,j);
        if(c<0){
          s << " x ";
        }
        else{
          s << " " << mat.data[c].v << " ";
        }
      }
      else{//i<j
        CPPL_INT c =mat.number(i,j);
        if(c<0){
          s << "{x}";
        }
        else{
          s << "{" << mat.data[c].v << "}";
        }
      }
    }
    s << std::endl;
  }
  
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline void dssmatrix::write(const char* filename) const
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
}

//=============================================================================
inline void dssmatrix::read(const char* filename)
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
  if( id != "dssmatrix" && id != "#dssmatrix" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not dssmatrix." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  //////// n ////////
  size_t size;
  s >> n >> size;
  resize(n);
  data.reserve(size);
  
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
