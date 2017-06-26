//=============================================================================
//! Samll Real Double-precision Symmetric Matrix Class
template<CPPL_INT n> class dsymatrix_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  double array[(n*(n+1))/2]; //!< 1D array to store vector data
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dsymatrix_small();
  inline explicit dsymatrix_small(const dsymatrix&);
  inline ~dsymatrix_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline dgematrix_small<n,n> to_dgematrix_small() const;
  inline dsymatrix to_dsymatrix() const;
  
  //////// io ////////
  inline double& operator()(const CPPL_INT& i, const CPPL_INT& j);
  inline double  operator()(const CPPL_INT& i, const CPPL_INT& j) const;
  inline dsymatrix_small<n>& set(const CPPL_INT&, const CPPL_INT&, const double&);
  template<CPPL_INT _n> inline friend std::ostream& operator<<(std::ostream&, const dsymatrix_small<_n>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// misc ////////
  inline dsymatrix_small<n>& zero();
  inline dsymatrix_small<n>& identity();
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT N> inline dsymatrix_small<N>& operator= (const dsymatrix_small<N>&);
  //////// += ////////
  template<CPPL_INT N> inline friend dsymatrix_small<N>& operator+=(dsymatrix_small<N>&, const dsymatrix_small<N>&);
  //////// -= ////////
  template<CPPL_INT N> inline friend dsymatrix_small<N>& operator-=(dsymatrix_small<N>&, const dsymatrix_small<N>&);
  //////// *= ////////
  template<CPPL_INT N> inline friend dsymatrix_small<N>& operator*=(dsymatrix_small<N>&, const dsymatrix_small<N>&);
  template<CPPL_INT N> inline friend dsymatrix_small<N>& operator*=(dsymatrix_small<N>&, const             double&);
  //////// /= ////////
  template<CPPL_INT N> inline friend dsymatrix_small<N>& operator/=(dsymatrix_small<N>&, const             double&);
  //////// unary ////////
  template<CPPL_INT N> inline friend const dsymatrix_small<N>& operator+(const dsymatrix_small<N>&);
  template<CPPL_INT N> inline friend dsymatrix_small< N > operator-(const dsymatrix_small<N>&);
  //////// + ////////
  template<CPPL_INT N> inline friend dgematrix_small<N,N> operator+(const dsymatrix_small<N>&, const dgematrix_small<N,N>&);
  template<CPPL_INT N> inline friend dsymatrix_small< N > operator+(const dsymatrix_small<N>&, const dsymatrix_small< N >&);
  //////// - ////////
  template<CPPL_INT N> inline friend dgematrix_small<N,N> operator-(const dsymatrix_small<N>&, const dgematrix_small<N,N>&);
  template<CPPL_INT N> inline friend dsymatrix_small< N > operator-(const dsymatrix_small<N>&, const dsymatrix_small< N >&);
  //////// * ////////
  template<CPPL_INT N> inline friend dcovector_small< N > operator*(const dsymatrix_small<N>&, const dcovector_small< N >&);
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator*(const dsymatrix_small<M>&, const dgematrix_small<M,N>&);
  template<CPPL_INT N> inline friend dgematrix_small<N,N> operator*(const dsymatrix_small<N>&, const dsymatrix_small< N >&);
  template<CPPL_INT N> inline friend dsymatrix_small< N > operator*(const dsymatrix_small<N>&, const               double&);
  //////// / ////////
  template<CPPL_INT N> inline friend dsymatrix_small< N > operator/(const dsymatrix_small<N>&, const               double&);
  //////// double ////////
  template<CPPL_INT N> inline friend dsymatrix_small< N > operator*(const             double&, const dsymatrix_small< N >&);
  //////// hadamerd ////////
  template<CPPL_INT N> inline friend dsymatrix_small< N >  hadamerd(const dsymatrix_small<N>&, const dgematrix_small<N,N>&);
  template<CPPL_INT N> inline friend dsymatrix_small< N >  hadamerd(const dsymatrix_small<N>&, const dsymatrix_small< N >&);
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline double det(const dsymat2&);
inline dsymat2 inv(const dsymat2&);
inline dsymat2 rotate(const dsymat2&, const double&);

inline double det(const dsymat3&);
inline dsymat3 inv(const dsymat3&);
inline dsymat3 rotate(const dsymat3&, const dquater&);
