//=============================================================================
/*! convert zgematrix_small to zgematrix */
template<CPPL_INT m, CPPL_INT n>
inline _zgematrix zgematrix_small<m,n>::to_zgematrix() const
{CPPL_VERBOSE_REPORT;
  zgematrix mat(m,n);
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      mat(i,j) =(*this)(i,j);
    }
  }
  return _(mat);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! operator() */
template<CPPL_INT m, CPPL_INT n>
inline comple& zgematrix_small<m,n>::operator()(const CPPL_INT& i, const CPPL_INT& j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  return array[i+m*j];
}

//=============================================================================
/*! operator() for const */
template<CPPL_INT m, CPPL_INT n>
inline comple zgematrix_small<m,n>::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i<0 || j<0 || m<=i || n<=j ){
    ERROR_REPORT;
    std::cerr << "The required component is out of the matrix size." << std::endl
              << "Your input is (" << i << "," << j << "), whereas the matrix size is " << m << "x" << n << "." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  return array[i+m*j];
}

//=============================================================================
/*! set */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& zgematrix_small<m,n>::set(const CPPL_INT& i, const CPPL_INT& j, const comple& v)
{CPPL_VERBOSE_REPORT;
  (*this)(i,j) =v;
  return *this;
}

//=============================================================================
/*! operator<< */
template<CPPL_INT m, CPPL_INT n>
inline std::ostream& operator<<(std::ostream& s, const zgematrix_small<m,n>& A)
{CPPL_VERBOSE_REPORT;
  s << std::setiosflags(std::ios::showpos);
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      s << " " << A(i,j);
    }
    s << std::endl;
  }
  return s;
}

//=============================================================================
/*! write to file */
template<CPPL_INT m, CPPL_INT n>
inline void zgematrix_small<m,n>::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  ofs << "#zgematrix" << " " << m << " " << n << std::endl;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      ofs << (*this)(i,j) << " ";
    }
    ofs << std::endl;
  }
  ofs.close();
}

//=============================================================================
/*! read from file */
template<CPPL_INT m, CPPL_INT n>
inline void zgematrix_small<m,n>::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s( filename );
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }
  
  std::string id;
  s >> id;
  if( id != "zgematrix" && id != "#zgematrix" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not zgematrix." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  CPPL_INT _m, _n;
  s >> _m >> _n;
  if(m!=_m || n!=_n){
    ERROR_REPORT;
    std::cerr << "Matrix size is invalid." << std::endl;
    exit(1);
  }
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++ ){
      s >> operator()(i,j);
    }
  }
  if(s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << "\"." << std::endl
              << "Most likely, there is not enough data components, or a linefeed code or space code is missing at the end of the last line." << std::endl;
    exit(1);
  }
  
  s >> id;
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
/*! return transposed zgematrix_small */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<n,m> t(const zgematrix_small<m,n>& A)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,m> X;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      X(j,i) =A(i,j);
    }
  }
  return X;
}

//=============================================================================
/*! return its trace */
template<CPPL_INT m, CPPL_INT n>
inline comple trace(const zgematrix_small<m,n>& A)
{CPPL_VERBOSE_REPORT;
  comple trace =comple(0.,0);
  
  const CPPL_INT imax =std::min(m,n);
  for(CPPL_INT i=0; i<imax; i++){
    trace +=A(i,i);
  }
  
  return trace;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zero */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& zgematrix_small<m,n>::zero()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<m*n; k++){
    array[k] =comple(0.,0.);
  }
  
  return *this;
}

//=============================================================================
/*! identity */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& zgematrix_small<m,n>::identity()
{CPPL_VERBOSE_REPORT;
  zero();
  
  const CPPL_INT kmax =std::min(m,n);
  for(CPPL_INT k=0; k<kmax; k++){
    (*this)(k,k) =1.;
  }
  
  return *this;
}

//=============================================================================
/*! return the j-th column vector */
template<CPPL_INT m, CPPL_INT n>
inline zcovector_small<m> zgematrix_small<m,n>::col(const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
  zcovector_small<m> vec;
  for(CPPL_INT i=0; i<m; i++){
    vec(i) =(*this)(i,j);
  }
  return vec;
}

//=============================================================================
/*! return the i-th row vector */
template<CPPL_INT m, CPPL_INT n>
inline zrovector_small<n> zgematrix_small<m,n>::row(const CPPL_INT& i) const
{CPPL_VERBOSE_REPORT;
  zrovector_small<n> vec;
  for(CPPL_INT j=0; j<n; j++){
    vec(j)=(*this)(i,j);
  }
  return vec;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix_small+=zgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& operator+=(zgematrix_small<m,n>& A, const zgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<m*n; k++){
    A.array[k] +=B.array[k];
  }
  return A;
}

//=============================================================================
/*! zgematrix_small-=zgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& operator-=(zgematrix_small<m,n>& A, const zgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<m*n; k++){
    A.array[k] -=B.array[k];
  }
  return A;
}

//=============================================================================
/*! zgematrix_small*=zgematrix_small operator */
template<CPPL_INT m, CPPL_INT l, CPPL_INT n>
inline zgematrix_small<m,n>& operator*=(zgematrix_small<m,l>& A, const zgematrix_small<l,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> X;
  X.zero();
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      for(CPPL_INT k=0; k<l; k++){
        X(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  A =X;
  return A;
}

//=============================================================================
/*! zgematrix_small*=double operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& operator*=(zgematrix_small<m,n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<m*n; k++){
    A.array[k] *=v;
  }
  return A;
}

//=============================================================================
/*! zgematrix_small*=comple operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& operator*=(zgematrix_small<m,n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<m*n; k++){
    A.array[k] *=v;
  }
  return A;
}

//=============================================================================
/*! zgematrix_small/=double operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& operator/=(zgematrix_small<m,n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<m*n; k++){
    A.array[k] /=v;
  }
  return A;
}
//=============================================================================
/*! zgematrix_small/=comple operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n>& operator/=(zgematrix_small<m,n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<m*n; k++){
    A.array[k] /=v;
  }
  return A;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! unary + operator */
template<CPPL_INT m, CPPL_INT n>
inline const zgematrix_small<m,n>& operator+(const zgematrix_small<m,n>& A)
{CPPL_VERBOSE_REPORT;
  return A;
}

//=============================================================================
/*! unary - operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator-(const zgematrix_small<m,n>& A)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> X;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      X(i,j) =-A(i,j);
    }
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix_small+zgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator+(const zgematrix_small<m,n>& A, const zgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      C(i,j) =A(i,j)+B(i,j);
    }
  }
  return C;
}

//=============================================================================
/*! zgematrix_small+zhematrix_small operator */
template<CPPL_INT n>
inline zgematrix_small<n,n> operator+(const zgematrix_small<n,n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> X;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<i; j++){
      X(i,j) =A(i,j)+B(i,j);
    }
    for(CPPL_INT j=i; j<n; j++){
      X(i,j) =A(i,j)+B(j,i);
    }
  }
  return X;
}

//=============================================================================
/*! zgematrix_small-zgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator-(const zgematrix_small<m,n>& A, const zgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      C(i,j)=A(i,j)-B(i,j);
    }
  }
  return C;
}

//=============================================================================
/*! zgematrix_small-zhematrix_small operator */
template<CPPL_INT n>
inline zgematrix_small<n,n> operator-(const zgematrix_small<n,n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> X;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      X(i,j)=A(i,j)-B(i,j);
    }
    for(CPPL_INT j=i+1; j<n; j++){
      X(i,j)=A(i,j)-B(j,i);
    }
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix_small*zcovector_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zcovector_small<m> operator*(const zgematrix_small<m,n>& A, const zcovector_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zcovector_small<m> C;
  C.zero();
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i) +=A(i,j)*B(j);
    }
  }
  return C;
}

//=============================================================================
/*! zgematrix_small*zgematrix_small operator */
template<CPPL_INT m, CPPL_INT l, CPPL_INT n>
inline zgematrix_small<m,n> operator*(const zgematrix_small<m,l>& A, const zgematrix_small<l,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  C.zero();
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      for(int k=0; k<l; k++){
        C(i,j) +=A(i,k)*B(k,j);
      }
    }
  }
  return C;
}

//=============================================================================
/*! zgematrix_small+zhematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator*(const zgematrix_small<m,n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> X;
  X.zero();
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      for(CPPL_INT k=0; k<j; k++){
        X(i,j) +=A(i,k)*B(j,k);
      }
      for(CPPL_INT k=j; k<n; k++){
        X(i,j) +=A(i,k)*B(k,j);
      }
    }
  }
  return X;
}

//=============================================================================
/*! zgematrix_small*double operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator*(const zgematrix_small<m,n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i,j) =A(i,j)*v;
    }
  }
  return C;
}

//=============================================================================
/*! zgematrix_small*comple operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator*(const zgematrix_small<m,n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i,j) =A(i,j)*v;
    }
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zgematrix_small/double operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator/(const zgematrix_small<m,n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i,j) =A(i,j)/v;
    }
  }
  return C;
}

//=============================================================================
/*! zgematrix_small/comple operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator/(const zgematrix_small<m,n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i,j) =A(i,j)/v;
    }
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! Hadamerd product */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> hadamerd(const zgematrix_small<m,n>& A, const zgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> C;
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      C(i,j) =A(i,j)*B(i,j);
    }
  }
  return C;
}

//=============================================================================
/*! Hadamerd product */
template<CPPL_INT n>
inline zgematrix_small<n,n> hadamerd(const zgematrix_small<n,n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> C;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      C(i,j) =A(i,j)*B(i,j);
    }
    for(CPPL_INT j=i+1; j<n; j++){
      C(i,j) =A(i,j)*conj(B(j,i));
    }
  }
  return C;
}
