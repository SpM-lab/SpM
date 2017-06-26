//=============================================================================
//! Samll Complex Double-precision General Dence Matrix Class
template<CPPL_INT m, CPPL_INT n> class zgematrix_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  comple array[m*n];
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline zgematrix_small();
  inline explicit zgematrix_small(const zgematrix&);
  inline ~zgematrix_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zgematrix to_zgematrix() const;
  
  //////// io ////////
  inline comple& operator()(const CPPL_INT& i, const CPPL_INT& j);
  inline comple  operator()(const CPPL_INT& i, const CPPL_INT& j) const;
  inline zgematrix_small<m,n>& set(const CPPL_INT& i, const CPPL_INT& j, const comple& v);
  template<CPPL_INT _m, CPPL_INT _n> inline friend std::ostream& operator<<(std::ostream&, const zgematrix_small<_m,_n>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// calc ////////
#ifndef _MSC_VER
  template<CPPL_INT _m, CPPL_INT _n> inline friend zgematrix_small<n,m> t(const zgematrix_small<m,n>&);
#endif//_MSC_VER
  
  //////// misc ////////
  inline zgematrix_small<m,n>& zero();
  inline zgematrix_small<m,n>& identity();
  inline zcovector_small<m> col(const CPPL_INT& j) const;
  inline zrovector_small<n> row(const CPPL_INT& i) const;

  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT M, CPPL_INT N> inline zgematrix_small<M,N>& operator= (const  zgematrix_small<M,N>&);
  //////// += ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N>& operator+=(zgematrix_small<M,N>&, const zgematrix_small<M,N>&);
  //////// -= ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N>& operator-=(zgematrix_small<M,N>&, const zgematrix_small<M,N>&);
  //////// *= ////////
  template<CPPL_INT M, CPPL_INT L, CPPL_INT N> inline friend zgematrix_small<M,N>& operator*=(zgematrix_small<M,L>&, const zgematrix_small<L,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N>& operator*=(zgematrix_small<M,N>&, const               double&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N>& operator*=(zgematrix_small<M,N>&, const               comple&);
  //////// /= ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N>& operator/=(zgematrix_small<M,N>&, const               double&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N>& operator/=(zgematrix_small<M,N>&, const               comple&);
  //////// unary ////////
  template<CPPL_INT M, CPPL_INT N> inline friend const zgematrix_small<M,N>& operator+(const zgematrix_small<M,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator-(const zgematrix_small<M,N>&);
  //////// + ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator+(const zgematrix_small<M,N>&, const zgematrix_small<M,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator+(const zgematrix_small<M,N>&, const zhematrix_small< N >&);
  //////// - ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator-(const zgematrix_small<M,N>&, const zgematrix_small<M,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator-(const zgematrix_small<M,N>&, const zhematrix_small< N >&);
  //////// * ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zcovector_small< M > operator*(const zgematrix_small<M,N>&, const zcovector_small< N >&);
  template<CPPL_INT M, CPPL_INT L, CPPL_INT N> inline friend zgematrix_small<M,N> operator*(const zgematrix_small<M,L>&, const zgematrix_small<L,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator*(const zgematrix_small<M,N>&, const zhematrix_small< N >&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator*(const zgematrix_small<M,N>&, const               double&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator*(const zgematrix_small<M,N>&, const               comple&);
  //////// / ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator/(const zgematrix_small<M,N>&, const               double&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator/(const zgematrix_small<M,N>&, const               comple&);
  //////// comple ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator*(const               double&, const zgematrix_small<M,N>&);
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator*(const               comple&, const zgematrix_small<M,N>&);
  //////// hadamerd ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N>  hadamerd(const zgematrix_small<M,N>&, const zgematrix_small<M,N>&);
  template<CPPL_INT N>         inline friend zgematrix_small<N,N>  hadamerd(const zgematrix_small<N,N>&, const zhematrix_small< N >&);
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline comple det(const zgemat2&);
inline zgemat2 inv(const zgemat2&);
inline comple det(const zgemat3&);
inline zgemat3 inv(const zgemat3&);
