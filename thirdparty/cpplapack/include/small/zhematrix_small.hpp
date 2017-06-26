//=============================================================================
//! Samll Complex Double-precision Symmetric Matrix Class
template<CPPL_INT n> class zhematrix_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  comple array[(n*(n+1))/2];
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline zhematrix_small();
  inline explicit zhematrix_small(const zhematrix&);
  inline ~zhematrix_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline zgematrix_small<n,n> to_zgematrix_small() const;
  inline zhematrix to_zhematrix() const;
  
  //////// io ////////
  inline comple& operator()(const CPPL_INT& i, const CPPL_INT& j);
  inline comple  operator()(const CPPL_INT& i, const CPPL_INT& j) const;
  inline zhematrix_small<n>& set(const CPPL_INT&, const CPPL_INT&, const comple&);
  template<CPPL_INT _n> inline friend std::ostream& operator<<(std::ostream&, const zhematrix_small<_n>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// misc ////////
  inline zhematrix_small<n>& zero();
  inline zhematrix_small<n>& identity();
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT N> inline zhematrix_small<N>& operator= (const zhematrix_small<N>&);
  //////// += ////////
  template<CPPL_INT N> inline friend zhematrix_small<N>& operator+=(zhematrix_small<N>&, const zhematrix_small<N>&);
  //////// -= ////////
  template<CPPL_INT N> inline friend zhematrix_small<N>& operator-=(zhematrix_small<N>&, const zhematrix_small<N>&);
  //////// *= ////////
  template<CPPL_INT N> inline friend zhematrix_small<N>& operator*=(zhematrix_small<N>&, const zhematrix_small<N>&);
  template<CPPL_INT N> inline friend zhematrix_small<N>& operator*=(zhematrix_small<N>&, const             double&);
  template<CPPL_INT N> inline friend zhematrix_small<N>& operator*=(zhematrix_small<N>&, const             comple&);
  //////// /= ////////
  template<CPPL_INT N> inline friend zhematrix_small<N>& operator/=(zhematrix_small<N>&, const             double&);
  template<CPPL_INT N> inline friend zhematrix_small<N>& operator/=(zhematrix_small<N>&, const             comple&);
  //////// unary ////////
  template<CPPL_INT N> inline friend const zhematrix_small<N>& operator+(const zhematrix_small<N>&);
  template<CPPL_INT N> inline friend zhematrix_small< N > operator-(const zhematrix_small<N>&);
  //////// + ////////
  template<CPPL_INT N> inline friend zgematrix_small<N,N> operator+(const zhematrix_small<N>&, const zgematrix_small<N,N>&);
  template<CPPL_INT N> inline friend zhematrix_small< N > operator+(const zhematrix_small<N>&, const zhematrix_small< N >&);
  //////// - ////////
  template<CPPL_INT N> inline friend zgematrix_small<N,N> operator-(const zhematrix_small<N>&, const zgematrix_small<N,N>&);
  template<CPPL_INT N> inline friend zhematrix_small< N > operator-(const zhematrix_small<N>&, const zhematrix_small< N >&);
  //////// * ////////
  template<CPPL_INT N> inline friend zcovector_small< N > operator*(const zhematrix_small<N>&, const zcovector_small< N >&);
  template<CPPL_INT N> inline friend zgematrix_small<N,N> operator*(const zhematrix_small<N>&, const zgematrix_small<N,N>&);
  template<CPPL_INT N> inline friend zgematrix_small<N,N> operator*(const zhematrix_small<N>&, const zhematrix_small< N >&);
  template<CPPL_INT N> inline friend zhematrix_small< N > operator*(const zhematrix_small<N>&, const               double&);
  template<CPPL_INT N> inline friend zhematrix_small< N > operator*(const zhematrix_small<N>&, const               comple&);
  //////// / ////////
  template<CPPL_INT N> inline friend zhematrix_small< N > operator/(const zhematrix_small<N>&, const               double&);
  template<CPPL_INT N> inline friend zhematrix_small< N > operator/(const zhematrix_small<N>&, const               comple&);
  //////// comple ////////
  template<CPPL_INT N> inline friend zhematrix_small< N > operator*(const             double&, const zhematrix_small< N >&);
  template<CPPL_INT N> inline friend zhematrix_small< N > operator*(const             comple&, const zhematrix_small< N >&);
  //////// hadamerd ////////
  template<CPPL_INT N> inline friend zgematrix_small<N,N>  hadamerd(const zhematrix_small<N>&, const zgematrix_small<N,N>&);
  template<CPPL_INT N> inline friend zhematrix_small< N >  hadamerd(const zhematrix_small<N>&, const zhematrix_small< N >&);  
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline comple det(const zhemat2&);
inline zhemat2 inv(const zhemat2&);
inline comple det(const zhemat3&);
inline zhemat3 inv(const zhemat3&);
