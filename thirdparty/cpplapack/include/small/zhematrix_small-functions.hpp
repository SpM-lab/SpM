//=============================================================================
/*! convert zhematrix_small to zgematrix_small */
template<CPPL_INT n>
inline zgematrix_small<n,n> zhematrix_small<n>::to_zgematrix_small() const
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> newmat;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0;   j<=i; j++){ newmat(i,j) =(*this)(i,j); }
    for(CPPL_INT j=i+1; j<n;  j++){ newmat(i,j) =(*this)(j,i); }
  }
  return newmat;
}

//=============================================================================
/*! convert zhematrix_small to zhematrix */
template<CPPL_INT n>
inline zhematrix zhematrix_small<n>::to_zhematrix() const
{CPPL_VERBOSE_REPORT;
  zhematrix newmat(n);
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      newmat(i,j) =(*this)(i,j);
    }
  }
  return newmat;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! operator() */
template<CPPL_INT n>
inline comple& zhematrix_small<n>::operator()(const CPPL_INT& i, const CPPL_INT& j)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(i<j){
    ERROR_REPORT;
    std::cerr << "i must be greater than or equal to j since dsymatrix_small is L strage. " << std::endl;
    std::cerr << "Your input was (" << i << ", " << j << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //const CPPL_INT I(max(i,j)), J(min(i,j)); return array[(I*(I+1))/2 +J];
  return array[(i*(i+1))/2 +j]; //l storage
}

//=============================================================================
/*! operator() for const */
template<CPPL_INT n>
inline comple zhematrix_small<n>::operator()(const CPPL_INT& i, const CPPL_INT& j) const
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(i<j){
    ERROR_REPORT;
    std::cerr << "i must be greater than or equal to j since dsymatrix_small is L strage. " << std::endl;
    std::cerr << "Your input was (" << i << ", " << j << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  //const CPPL_INT I(max(i,j)), J(min(i,j)); return array[(I*(I+1))/2 +J];
  return array[(i*(i+1))/2 +j];
}

//=============================================================================
/*! set */
template<CPPL_INT n>
inline zhematrix_small<n>& zhematrix_small<n>::set(const CPPL_INT& i, const CPPL_INT& j, const comple& v)
{CPPL_VERBOSE_REPORT;
  (*this)(i,j)=v;
  return *this;
}

//=============================================================================
/*! operator<< */
template<CPPL_INT n>
inline std::ostream& operator<<(std::ostream& s, const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  s << std::setiosflags(std::ios::showpos);
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0;   j<=i; j++){ s << " " << A(i,j) << " "<< std::flush; }
    for(CPPL_INT j=i+1; j<n;  j++){ s << "{" << A(j,i) << "}" << std::flush; }
    s << std::endl;
  }
  return s;
}

//=============================================================================
/*! write to file */
template<CPPL_INT n>
inline void zhematrix_small<n>::write(const char* filename) const
{CPPL_VERBOSE_REPORT;
  std::ofstream ofs(filename, std::ios::trunc);
  ofs.setf(std::cout.flags());
  ofs.precision(std::cout.precision());
  ofs.width(std::cout.width());
  ofs.fill(std::cout.fill());
  ofs << "#zhematrix" << " " << n << std::endl;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      ofs << (*this)(i,j) << " ";
    }
    ofs << std::endl;
  }
  ofs.close();
}

//=============================================================================
/*! read from file */
template<CPPL_INT n>
inline void zhematrix_small<n>::read(const char* filename)
{CPPL_VERBOSE_REPORT;
  std::ifstream s(filename);
  if(!s){
    ERROR_REPORT;
    std::cerr << "The file \"" << filename << "\" can not be opened." << std::endl;
    exit(1);
  }
  
  std::string id;
  s >> id;
  if( id != "zhematrix" && id != "#zhematrix" ){
    ERROR_REPORT;
    std::cerr << "The type name of the file \"" << filename << "\" is not zhematrix." << std::endl
              << "Its type name was " << id << " ." << std::endl;
    exit(1);
  }
  
  CPPL_INT _n;
  s >> _n;
  if(n!=_n){
    ERROR_REPORT;
    std::cerr << "Matrix size is invalid." << std::endl;
    exit(1);
  }
  
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++ ){
      s >> operator()(i,j);
    }
  }
  if(s.eof()){
    ERROR_REPORT;
    std::cerr << "There is something is wrong with the file \"" << filename << "\"." << std::endl
              << "Most likely, there is a lack of data components, or a linefeed code or space code is missing at the end of the last line." << std::endl;
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
/*! zero */
template<CPPL_INT n>
inline zhematrix_small<n>& zhematrix_small<n>::zero()
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      (*this)(i,j) =0.;
    }
  }
  return *this;
}

//=============================================================================
/*! identity */
template<CPPL_INT n>
inline zhematrix_small<n>& zhematrix_small<n>::identity()
{CPPL_VERBOSE_REPORT;
  zero();
  for(CPPL_INT k=0; k<n; k++){
    (*this)(k,k) =1.;
  }
  return *this;
}

//=============================================================================
/*! return its trace */
template<CPPL_INT n>
inline comple trace(const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  comple trace =comple(0.,0.);
  for(CPPL_INT i=0; i<n; i++){
    trace +=A(i,i);
  }
  return trace;
}

//=============================================================================
/*! find index of the maximum component */
template<CPPL_INT n>
inline void idamax(CPPL_INT& I, CPPL_INT& J, const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  double max(-1.);
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      if( max<fabs(A(i,j)) ){
        I=i;
        J=j;
        max =fabs(A(i,j));
      }
    }
  }
  return;  
}

//=============================================================================
/*! return the maximum component */
template<CPPL_INT n>
inline comple damax(const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  CPPL_INT i(0),j(0);
  idamax(i,j,A);
  return A(i,j);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix_small+=zhematrix_small operator */
template<CPPL_INT n>
inline zhematrix_small<n>& operator+=(zhematrix_small<n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<(n*(n+1))/2; k++){
    A.array[k] +=B.array[k];
  }
  return A;
}

//=============================================================================
/*! zhematrix_small-=zhematrix_small operator */
template<CPPL_INT n>
inline zhematrix_small<n>& operator-=(zhematrix_small<n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT k=0; k<(n*(n+1))/2; k++){
    A.array[k] -=B.array[k];
  }
  return A;
}

//=============================================================================
/*! zhematrix_small*=double operator */
template<CPPL_INT n>
inline zhematrix_small<n>& operator*=(zhematrix_small<n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      A(i,j) *=v;
    }
  }
  return A;
}

//=============================================================================
/*! zhematrix_small*=comple operator */
template<CPPL_INT n>
inline zhematrix_small<n>& operator*=(zhematrix_small<n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      A(i,j) *=v;
    }
  }
  return A;
}

//=============================================================================
/*! zhematrix_small/=double operator */
template<CPPL_INT n>
inline zhematrix_small<n>& operator/=(zhematrix_small<n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      A(i,j) /=v;
    }
  }
  return A;
}

//=============================================================================
/*! zhematrix_small/=comple operator */
template<CPPL_INT n>
inline zhematrix_small<n>& operator/=(zhematrix_small<n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      A(i,j) /=v;
    }
  }
  return A;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! unary */
template<CPPL_INT n>
inline const zhematrix_small<n>& operator+(const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  return A;
}

//=============================================================================
/*! unary */
template<CPPL_INT n>
inline zhematrix_small<n> operator-(const zhematrix_small<n>& A)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> X;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      X(i,j) =-A(i,j);
    }
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix_small+zgematrix_small operator */
template<CPPL_INT n>
inline zgematrix_small<n,n> operator+(const zhematrix_small<n>& A, const zgematrix_small<n,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> X;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<i; j++){
      X(i,j) =A(i,j)+B(i,j);
    }
    for(CPPL_INT j=i; j<n; j++){
      X(i,j) =A(j,i)+B(i,j);
    }
  }
  return X;
}

//=============================================================================
/*! zhematrix_small+zhematrix_small operator */
template<CPPL_INT n>
inline zhematrix_small<n> operator+(const zhematrix_small<n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> X;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      X(i,j) =A(i,j)+B(i,j);
    }
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix_small-zgematrix_small operator */
template<CPPL_INT n>
inline zgematrix_small<n,n> operator-(const zhematrix_small<n>& A, const zgematrix_small<n,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> X;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<i; j++){
      X(i,j) =A(i,j)-B(i,j);
    }
    for(CPPL_INT j=i; j<n; j++){
      X(i,j) =A(j,i)-B(i,j);
    }
  }
  return X;
}

//=============================================================================
/*! zhematrix_small-zhematrix_small operator */
template<CPPL_INT n>
inline zhematrix_small<n> operator-(const zhematrix_small<n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> X;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      X(i,j) =A(i,j)-B(i,j);
    }
  }
  return X;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix_small*zcovector_small operator */
template<CPPL_INT n>
inline zcovector_small<n> operator*(const zhematrix_small<n>& A, const zcovector_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zcovector_small<n> C;
  C.zero();
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<i; j++){
      C(i) +=A(i,j)*B(j);
    }
    for(CPPL_INT j=i; j<n; j++){
      C(i) +=A(j,i)*B(j);
    }
  }
  return C;
}

//=============================================================================
/*! zhematrix_small*zgematrix_small operator */
template<CPPL_INT m, CPPL_INT n>
inline zgematrix_small<m,n> operator*(const zhematrix_small<m>& A, const zgematrix_small<m,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<m,n> X;
  X.zero();
  for(CPPL_INT i=0; i<m; i++){
    for(CPPL_INT j=0; j<n; j++){
      for(CPPL_INT k=0; k<i; k++){
        X(i,j) +=A(i,k)*B(k,j);
      }
      for(CPPL_INT k=i; k<m; k++){
        X(i,j) +=A(k,i)*B(k,j);
      }
    }
  }
  return X;
}

//=============================================================================
/*! zhematrix_small*zhematrix_small operator */
template<CPPL_INT n>
inline zgematrix_small<n,n> operator*(const zhematrix_small<n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> X;
  X.zero();
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<i; j++){
      for(CPPL_INT k=0; k<j; k++){
        X(i,j) +=A(i,k)*B(j,k);
      }
      for(CPPL_INT k=j; k<i; k++){
        X(i,j) +=A(i,k)*B(k,j);
      }
      for(CPPL_INT k=i; k<n; k++){
        X(i,j) +=A(k,i)*B(k,j);
      }
    }
    for(CPPL_INT j=i; j<n; j++){
      for(CPPL_INT k=0; k<i; k++){
        X(i,j) +=A(i,k)*B(j,k);
      }
      for(CPPL_INT k=i; k<j; k++){
        X(i,j) +=A(k,i)*B(j,k);
      }
      for(CPPL_INT k=j; k<n; k++){
        X(i,j) +=A(k,i)*B(k,j);
      }
    }
  }
  return X;
}

//=============================================================================
/*! zhematrix_small*double operator */
template<CPPL_INT n>
inline zhematrix_small<n> operator*(const zhematrix_small<n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> C;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      C(i,j) =A(i,j)*v;
    }
  }
  return C;
}

//=============================================================================
/*! zhematrix_small*comple operator */
template<CPPL_INT n>
inline zhematrix_small<n> operator*(const zhematrix_small<n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> C;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      C(i,j) =A(i,j)*v;
    }
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! zhematrix_small/double operator */
template<CPPL_INT n>
inline zhematrix_small<n> operator/(const zhematrix_small<n>& A, const double& v)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> C;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      C(i,j) =A(i,j)/v;
    }
  }
  return C;
}

//=============================================================================
/*! zhematrix_small/comple operator */
template<CPPL_INT n>
inline zhematrix_small<n> operator/(const zhematrix_small<n>& A, const comple& v)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> C;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      C(i,j) =A(i,j)/v;
    }
  }
  return C;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! Hadamard product */
template<CPPL_INT n>
inline zgematrix_small<n,n> hadamard(const zhematrix_small<n>& A, const zgematrix_small<n,n>& B)
{CPPL_VERBOSE_REPORT;
  zgematrix_small<n,n> C;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      C(i,j) =A(i,j)*B(i,j);
    }
    for(CPPL_INT j=i+1; j<n; j++){
      C(i,j) =conj(A(j,i))*B(i,j);
    }
  }
  return C;
}

//=============================================================================
/*! Hadamard product */
template<CPPL_INT n>
inline zhematrix_small<n> hadamard(const zhematrix_small<n>& A, const zhematrix_small<n>& B)
{CPPL_VERBOSE_REPORT;
  zhematrix_small<n> C;
  for(CPPL_INT i=0; i<n; i++){
    for(CPPL_INT j=0; j<=i; j++){
      C(i,j) =A(i,j)*B(i,j);
    }
  }
  return C;
}
