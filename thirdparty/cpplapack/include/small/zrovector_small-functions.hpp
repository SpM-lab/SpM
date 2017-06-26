//=============================================================================
/*! convert zrovector_small to zrovector */
template<CPPL_INT l>
inline _zrovector zrovector_small<l>::to_zrovector() const
{CPPL_VERBOSE_REPORT;
  zrovector vec(l);
  for(CPPL_INT k=0; k<l; k++){
    vec(k) =(*this)(k);
  }
  return _(vec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! operator() */
template<CPPL_INT l>
inline comple& zrovector_small<l>::operator()(const CPPL_INT& k)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( k<0 || l<=k ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the vector size." << std::endl
              << "Your input is (" << k << "), whereas the vector size is " << l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  return array[k];
}

//=============================================================================
/*! operator() for const */
template<CPPL_INT l>
inline comple zrovector_small<l>::operator()(const CPPL_INT& k) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( k<0 || l<=k ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the vector size." << std::endl
              << "Your input is (" << k << "), whereas the vector size is " << l << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  return array[k];
}

//=============================================================================
/*! set */
template<CPPL_INT l>
inline zrovector_small<l>& zrovector_small<l>::set(const CPPL_INT& k, const comple& v)
{CPPL_VERBOSE_REPORT;
  (*this)(k) =v;
  return *this;
}

//=============================================================================
/*! operator<< */
template<CPPL_INT l>
inline std::ostream& operator<<(std::ostream& s, const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  s << std::setiosflags(std::ios::showpos);
  for(CPPL_INT i=0; i<l; i++){
    s << " " << A(i) << std::flush;
  }
  s << std::endl;
  return s;
}

//=============================================================================
/*! write to file */
template<CPPL_INT l>
inline void zrovector_small<l>::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#zrovector" << " " << l << std::endl;
  for(CPPL_INT k=0; k<l; k++){
    ofs << (*this)(k) << std::endl;
  }
  ofs.close();
}

//=============================================================================
/*! read from file */
template<CPPL_INT l>
inline void zrovector_small<l>::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }
  
  std::string id;
  s >> id;
  if( id != "zrovector" && id != "#zrovector" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not zrovector." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  CPPL_INT _l;
  s >> _l;
  if(l!=_l){
    ERROR_REPORT;
    std::cerr << "Matrix size is invalid." << std::endl;
    exit(1);
  }
  for(CPPL_INT k=0; k<l; k++){
    s >> (*this)(k);
  }
  if(s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << "\"." << std::endl
              << "Most likely, there is not enough data components, or a linefeed code or space code is missing at the end of the last line." << std::endl;
    exit(1);
  }
  
  s >> id;//tmp
  if(!s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << "\"." << std::endl
              << "Most likely, there are extra data components." << std::endl;
    exit(1);
  }
  
  s.close();    
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! return transposed zrovector_small */
template<CPPL_INT n>
inline zcovector_small<n> t(const zrovector_small<n>& A)
{CPPL_VERBOSE_REPORT;
  zcovector_small<n> X;
  for(CPPL_INT i=0; i<n; i++){
    X(i)=A(i);
  }
  return X;
}

//=============================================================================
/*! return its 2-norm */
template<CPPL_INT l>
inline double nrm2(const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  double v(0.);
  for(CPPL_INT i=0; i<l; i++){
    v+=A(i)*A(i);
  }
  return std::sqrt(v);
}

//=============================================================================
/*! find index of the maximum component */
template<CPPL_INT l>
inline void idamax(CPPL_INT& K, const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  double max(-1.);
  for(int k=0; k<l; k++){
    if( max<fabs(A(k)) ){
      K=k;
      max =fabs(A(k));
    }
  }
  return;
}

//=============================================================================
/*! return the maximum component */
template<CPPL_INT l>
inline comple damax(const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  CPPL_INT k(0);
  idamax(k,A);
  return A(k);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zero */
template<CPPL_INT l>
inline zrovector_small<l>& zrovector_small<l>::zero()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<l; k++){
    array[k] =0.;
  }
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector_small+=zrovector_small operator */
template<CPPL_INT l>
inline zrovector_small<l>& operator+=(zrovector_small<l>& A, const zrovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) +=B(i);
  }
  return A;
}

//=============================================================================
/*! zrovector_small-=zrovector_small operator */
template<CPPL_INT l>
inline zrovector_small<l>& operator-=(zrovector_small<l>& A, const zrovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) -=B(i);
  }
  return A;
}

//=============================================================================
/*! zrovector_small*=double operator */
template<CPPL_INT l>
inline zrovector_small<l>& operator*=(zrovector_small<l>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) *=v;
  }
  return A;
}

//=============================================================================
/*! zrovector_small*=comple operator */
template<CPPL_INT l>
inline zrovector_small<l>& operator*=(zrovector_small<l>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) *=v;
  }
  return A;
}

//=============================================================================
/*! zrovector_small/=double operator */
template<CPPL_INT l>
inline zrovector_small<l>& operator/=(zrovector_small<l>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) /=v;
  }
  return A;
}

//=============================================================================
/*! zrovector_small/=comple operator */
template<CPPL_INT l>
inline zrovector_small<l>& operator/=(zrovector_small<l>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) /=v;
  }
  return A;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! unary + operator */
template<CPPL_INT l>
inline const zrovector_small<l>& operator+(const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  return A;
}

//=============================================================================
/*! unary - operator */
template<CPPL_INT l>
inline zrovector_small<l> operator-(const zrovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =-A(i);
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector_small+zrovector_small operator */
template<CPPL_INT l>
inline zrovector_small<l> operator+(const zrovector_small<l>& A, const zrovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =A(i)+B(i);
  }
  return X;
}

//=============================================================================
/*! zrovector_small-zrovector_small operator */
template<CPPL_INT l>
inline zrovector_small<l> operator-(const zrovector_small<l>& A, const zrovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =A(i)-B(i);
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector_small*zcovector_small operator */
template<CPPL_INT l>
inline comple operator*(const zrovector_small<l>& A, const zcovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  comple x =0.;
  for(CPPL_INT i=0; i<l; i++){
    x +=A(i)*B(i);
  }
  return x;
}

//=============================================================================
/*! zrovector_small*zgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zrovector_small<n> operator*(const zrovector_small<m>& A, const zgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  zrovector_small<n> C;
  C.zero();
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=0; i<m; i++){
      C(j) +=A(i)*B(i,j);
    }
  }
  return C;
}

//=============================================================================
/*! zrovector_small*zhematrix_small operator */
template<CPPL_INT l>
inline zrovector_small<l> operator*(const zrovector_small<l>& A, const zhematrix_small<l>& B)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> C;
  C.zero();
  for(CPPL_INT j=0; j<l; j++){
    for(CPPL_INT i=0; i<j; i++){
      C(j) +=A(i)*B(j,i);
    }
    for(CPPL_INT i=j; i<l; i++){
      C(j) +=A(i)*B(i,j);
    }    
  }
  return C;
}

//=============================================================================
/*! zrovector_small*double operator */
template<CPPL_INT l>
inline zrovector_small<l> operator*(const zrovector_small<l>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)*v;
  }
  return C;
}

//=============================================================================
/*! zrovector_small*comple operator */
template<CPPL_INT l>
inline zrovector_small<l> operator*(const zrovector_small<l>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)*v;
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zrovector_small/double operator */
template<CPPL_INT l>
inline zrovector_small<l> operator/(const zrovector_small<l>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)/v;
  }
  return C;
}

//=============================================================================
/*! zrovector_small/comple operator */
template<CPPL_INT l>
inline zrovector_small<l> operator/(const zrovector_small<l>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)/v;
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! Hadamard product */
template<CPPL_INT l>
inline zrovector_small<l> hadamard(const zrovector_small<l>& A, const zrovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  zrovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)*B(i);
  }
  return C;
}
