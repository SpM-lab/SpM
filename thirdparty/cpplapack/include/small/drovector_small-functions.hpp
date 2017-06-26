//=============================================================================
/*! convert drovector_small to drovector */
template<CPPL_INT l>
inline _drovector drovector_small<l>::to_drovector() const
{CPPL_VERBOSE_REPORT;
  drovector vec(l);
  for(CPPL_INT k=0; k<l; k++){
    vec(k)=(*this)(k);
  }
  return _(vec);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! operator() */
template<CPPL_INT l>
inline double& drovector_small<l>::operator()(const CPPL_INT& k)
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
inline double drovector_small<l>::operator()(const CPPL_INT& k) const
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
/*! set function */
template<CPPL_INT l>
inline drovector_small<l>& drovector_small<l>::set(const CPPL_INT& k, const double& v)
{CPPL_VERBOSE_REPORT;
  (*this)(k) =v;
  return *this;
}

//=============================================================================
/*! operator<< */
template<CPPL_INT l>
inline std::ostream& operator<<(std::ostream& s, const drovector_small<l>& A)
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
inline void drovector_small<l>::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#drovector" << " " << l << std::endl;
  for(CPPL_INT k=0; k<l; k++){
    ofs << (*this)(k) << " ";
  }
  ofs << std::endl;
  ofs.close();
}

//=============================================================================
/*! read from file */
template<CPPL_INT l>
inline void drovector_small<l>::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }
  
  std::string id;
  s >> id;
  if( id != "drovector" && id != "#drovector" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not drovector." << std::endl
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
/*! return transposed column vector  */
template<CPPL_INT n>
inline dcovector_small<n> t(const drovector_small<n>& A)
{CPPL_VERBOSE_REPORT;
  dcovector_small<n> X;
  for(CPPL_INT i=0; i<n; i++){
    X(i)=A(i);
  }
  return X;
}

//=============================================================================
/*!  */
template<CPPL_INT l>
inline double nrm2(const drovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  double v(0.);
  for(CPPL_INT i=0; i<l; i++){
    v+=A(i)*A(i);
  }
  return std::sqrt(v);
}

//=============================================================================
/*!  */
template<CPPL_INT l>
inline void idamax(CPPL_INT& K, const drovector_small<l>& A)
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
/*!  */
template<CPPL_INT l>
inline double damax(const drovector_small<l>& A)
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
/*!  */
template<CPPL_INT l>
inline drovector_small<l>& drovector_small<l>::zero()
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
/*!  */
template<CPPL_INT l>
inline drovector_small<l>& operator+=(drovector_small<l>& A, const drovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) +=B(i);
  }
  return A;
}

//=============================================================================
/*!  */
template<CPPL_INT l>
inline drovector_small<l>& operator-=(drovector_small<l>& A, const drovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) -=B(i);
  }
  return A;
}

//=============================================================================
/*!  */
template<CPPL_INT l>
inline drovector_small<l>& operator*=(drovector_small<l>& A, const double& d)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) *=d;
  }
  return A;
}

//=============================================================================
/*!  */
template<CPPL_INT l>
inline drovector_small<l>& operator/=(drovector_small<l>& A, const double& d)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) /=d;
  }
  return A;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! unary */
template<CPPL_INT l>
inline const drovector_small<l>& operator+(const drovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  return A;
}

//=============================================================================
/*! unary */
template<CPPL_INT l>
inline drovector_small<l> operator-(const drovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =-A(i);
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*!  */
template<CPPL_INT l>
inline drovector_small<l> operator+(const drovector_small<l>& A, const drovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =A(i)+B(i);
  }
  return X;
}

//=============================================================================
/*!  */
template<CPPL_INT l>
inline drovector_small<l> operator-(const drovector_small<l>& A, const drovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =A(i)-B(i);
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*!  */
template<CPPL_INT l>
inline double operator*(const drovector_small<l>& A, const dcovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  double x =0.;
  for(CPPL_INT i=0; i<l; i++){
    x +=A(i)*B(i);
  }
  return x;
}

//=============================================================================
/*!  */
template<CPPL_INT m, CPPL_INT n>
inline drovector_small<n> operator*(const drovector_small<m>& A, const dgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  drovector_small<n> C;
  C.zero();
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=0; i<m; i++){
      C(j) +=A(i)*B(i,j);
    }
  }
  return C;
}

//=============================================================================
/*!  */
template<CPPL_INT l>
inline drovector_small<l> operator*(const drovector_small<l>& A, const dsymatrix_small<l>& B)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> C;
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
/*!  */
template<CPPL_INT l>
inline drovector_small<l> operator*(const drovector_small<l>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)*v;
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*!  */
template<CPPL_INT l>
inline drovector_small<l> operator/(const drovector_small<l>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)/v;
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! drovector_small%drovector_small (inner product) operator */
template<CPPL_INT l>
inline double operator%(const drovector_small<l>& A, const drovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  double v(0.);
  for(CPPL_INT i=0; i<l; i++){
    v +=A(i)*B(i);
  }
  return v;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! Hadamard product */
template<CPPL_INT l>
inline drovector_small<l> hadamard(const drovector_small<l>& A, const drovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)*B(i);
  }
  return C;
}
