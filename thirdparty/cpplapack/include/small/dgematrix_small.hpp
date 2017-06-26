//=============================================================================
//! Samll Real Double-precision General Dence Matrix Class
template<CPPL_INT m, CPPL_INT n> class dgematrix_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  double array[m*n]; //!< 1D array to store vector data
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dgematrix_small();
  inline explicit dgematrix_small(const dgematrix&);
  inline ~dgematrix_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _dgematrix to_dgematrix() const;
  
  //////// io ////////
  inline double& operator()(const CPPL_INT& i, const CPPL_INT& j);
  inline double  operator()(const CPPL_INT& i, const CPPL_INT& j) const;
  inline dgematrix_small<m,n>& set(const CPPL_INT& i, const CPPL_INT& j, const double& v);
  template<CPPL_INT _m, CPPL_INT _n> inline friend std::ostream& operator<<(std::ostream&, const dgematrix_small<_m,_n>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// calc ////////
#ifndef _MSC_VER
  template<CPPL_INT _m, CPPL_INT _n> inline friend dgematrix_small<n,m> t(const dgematrix_small<m,n>&);
  template<CPPL_INT _m, CPPL_INT _n> inline friend void idamax(CPPL_INT&, CPPL_INT&, const dgematrix_small&);
  template<CPPL_INT _m, CPPL_INT _n> inline friend double damax(const dgematrix_small&);
#endif//_MSC_VER
  
  //////// misc ////////
  inline dgematrix_small<m,n>& zero();
  inline dgematrix_small<m,n>& identity();
  inline dcovector_small<m> col(const CPPL_INT& j) const;
  inline drovector_small<n> row(const CPPL_INT& i) const;

  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT M, CPPL_INT N> inline dgematrix_small<M,N>& operator=(const dgematrix_small<M,N>&);
  //////// += ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N>& operator+=(dgematrix_small<M,N>&, const dgematrix_small<M,N>&);
  //////// -= ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N>& operator-=(dgematrix_small<M,N>&, const dgematrix_small<M,N>&);
  //////// *= ////////
  template<CPPL_INT M, CPPL_INT L, CPPL_INT N> inline friend dgematrix_small<M,N>& operator*=(dgematrix_small<M,L>&, const dgematrix_small<L,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N>& operator*=(dgematrix_small<M,N>&, const               double&);
  //////// /= ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N>& operator/=(dgematrix_small<M,N>&, const               double&);
  //////// unary ////////
  template<CPPL_INT M, CPPL_INT N> inline friend const dgematrix_small<M,N>& operator+(const dgematrix_small<M,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator-(const dgematrix_small<M,N>&);
  //////// + ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator+(const dgematrix_small<M,N>&, const dgematrix_small<M,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator+(const dgematrix_small<M,N>&, const dsymatrix_small< N >&);
  //////// - ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator-(const dgematrix_small<M,N>&, const dgematrix_small<M,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator-(const dgematrix_small<M,N>&, const dsymatrix_small< N >&);
  //////// * ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dcovector_small< M > operator*(const dgematrix_small<M,N>&, const dcovector_small< N >&);
  template<CPPL_INT M, CPPL_INT L, CPPL_INT N> inline friend dgematrix_small<M,N> operator*(const dgematrix_small<M,L>&, const dgematrix_small<L,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator*(const dgematrix_small<M,N>&, const dsymatrix_small< N >&);
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator*(const dgematrix_small<M,N>&, const               double&);
  //////// / ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator/(const dgematrix_small<M,N>&, const               double&);
  //////// double ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator*(const               double&, const dgematrix_small<M,N>&);
  //////// hadamerd ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N>  hadamerd(const dgematrix_small<M,N>&, const dgematrix_small<M,N>&);
  template<CPPL_INT N>         inline friend dgematrix_small<N,N>  hadamerd(const dgematrix_small<N,N>&, const dsymatrix_small< N >&);
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline double det(const dgemat2&);
inline dgemat2 inv(const dgemat2&);
inline dgemat2 rotate(const dgemat2&, const double&);
inline dgemat2 t2m(const double&);

inline double det(const dgemat3&);
inline dgemat3 inv(const dgemat3&);
inline dgemat3 rotate(const dgemat3&, const dquater&);
inline dquater m2q(const dgemat3&);

inline double det(const dgemat4&);
inline dgemat4 inv(const dgemat4&);
