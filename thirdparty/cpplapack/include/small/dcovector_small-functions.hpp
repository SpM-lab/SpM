//=============================================================================
/*! convert dcovector_small to dcovector */
template<CPPL_INT l>
inline _dcovector dcovector_small<l>::to_dcovector() const
{CPPL_VERBOSE_REPORT;
  dcovector vec(l);
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
inline double& dcovector_small<l>::operator()(const CPPL_INT& k)
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
inline double dcovector_small<l>::operator()(const CPPL_INT& k) const
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
inline dcovector_small<l>& dcovector_small<l>::set(const CPPL_INT& k, const double& v)
{CPPL_VERBOSE_REPORT;
  (*this)(k) =v;
  return *this;
}

//=============================================================================
/*! operator<< */
template<CPPL_INT l>
inline std::ostream& operator<<(std::ostream& s, const dcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  s << std::setiosflags(std::ios::showpos);
  for(CPPL_INT i=0; i<l; i++){
    s << A(i) << std::endl;
  }
  return s;
}

//=============================================================================
/*! write to file */
template<CPPL_INT l>
inline void dcovector_small<l>::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  
  ofs << "#dcovector" << " " << l << std::endl;
  for(CPPL_INT k=0; k<l; k++){
    ofs << (*this)(k) << std::endl;
  }
  ofs.close();
}

//=============================================================================
/*! read from file */
template<CPPL_INT l>
inline void dcovector_small<l>::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }
  
  std::string id;
  s >> id;
  if( id != "dcovector" && id != "#dcovector" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not dcovector." << std::endl
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
/*! return transposed drovector_small */
template<CPPL_INT l>
inline drovector_small<l> t(const dcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  drovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =A(i);
  }
  return X;
}

//=============================================================================
/*! return its 2-norm */
template<CPPL_INT l>
inline double nrm2(const dcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  double v(0);
  for(CPPL_INT i=0; i<l; i++){
    v+=A(i)*A(i);
  }
  return std::sqrt(v);
}

//=============================================================================
/*! return index of the maximum component */
template<CPPL_INT l>
inline void idamax(CPPL_INT& K, const dcovector_small<l>& A)
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
inline double damax(const dcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  CPPL_INT k(0);
  idamax(k,A);
  return A(k);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zero */
template<CPPL_INT l>
inline dcovector_small<l>& dcovector_small<l>::zero()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<l; k++){
    array[k] =0.;
  }
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector_small+=dcovector_small operator */
template<CPPL_INT l>
inline dcovector_small<l>& operator+=(dcovector_small<l>& A, const dcovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) +=B(i);
  }
  return A;
}

//=============================================================================
/*! dcovector_small-=dcovector_small operator */
template<CPPL_INT l>
inline dcovector_small<l>& operator-=(dcovector_small<l>& A, const dcovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) -=B(i);
  }
  return A;
}

//=============================================================================
/*! dcovector_small*=double operator */
template<CPPL_INT l>
inline dcovector_small<l>& operator*=(dcovector_small<l>& A, const double& d)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<l; i++){
    A(i) *=d;
  }
  return A;
}

//=============================================================================
/*! dcovector_small/=double operator */
template<CPPL_INT l>
inline dcovector_small<l>& operator/=(dcovector_small<l>& A, const double& d)
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
/*! unary + operator*/
template<CPPL_INT l>
inline const dcovector_small<l>& operator+(const dcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  return A;
}

//=============================================================================
/*! unary - operator*/
template<CPPL_INT l>
inline dcovector_small<l> operator-(const dcovector_small<l>& A)
{CPPL_VERBOSE_REPORT;
  dcovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =-A(i);
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector_small+dcovector_small operator */
template<CPPL_INT l>
inline dcovector_small<l> operator+(const dcovector_small<l>& A, const dcovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  dcovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =A(i)+B(i);
  }
  return X;
}

//=============================================================================
/*! dcovector_small-dcovector_small operator */
template<CPPL_INT l>
inline dcovector_small<l> operator-(const dcovector_small<l>& A, const dcovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  dcovector_small<l> X;
  for(CPPL_INT i=0; i<l; i++){
    X(i) =A(i)-B(i);
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector_small*double operator */
template<CPPL_INT n>
inline dcovector_small<n> operator*(const dcovector_small<n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  dcovector_small<n> C;
  for(CPPL_INT i=0; i<n; i++){
    C(i) =A(i)*v;
  }
  return C;
}

//=============================================================================
/*! dcovector_small*drovector_small operator */
template<CPPL_INT m, CPPL_INT n>
inline dgematrix_small<m,n> operator*(const dcovector_small<m>& A, const drovector_small<n>& B)
{CPPL_VERBOSE_REPORT;
  dgematrix_small<m,n> mat;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      mat(i,j) =A(i)*B(j);
    }
  }
  return mat;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector_small/double operator */
template<CPPL_INT n>
inline dcovector_small<n> operator/(const dcovector_small<n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  dcovector_small<n> C;
  for(CPPL_INT i=0; i<n; i++){
    C(i) =A(i)/v;
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! dcovector_small%dcovector_small (inner product) operator */
template<CPPL_INT l>
inline double operator%(const dcovector_small<l>& A, const dcovector_small<l>& B)
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
inline dcovector_small<l> hadamard(const dcovector_small<l>& A, const dcovector_small<l>& B)
{CPPL_VERBOSE_REPORT;
  dcovector_small<l> C;
  for(CPPL_INT i=0; i<l; i++){
    C(i) =A(i)*B(i);
  }
  return C;
}

