//=============================================================================
//! Samll Complex Double-precision Row Vector Class
template<CPPL_INT l> class zrovector_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  comple array[l];
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline zrovector_small();
  inline explicit zrovector_small(const zrovector&);
  inline zrovector_small(const comple&, const comple&);
  inline zrovector_small(const comple&, const comple&, const comple&);
  inline ~zrovector_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zrovector to_zrovector() const;
  
  //////// io ////////
  inline comple& operator()(const CPPL_INT&);
  inline comple  operator()(const CPPL_INT&) const;
  inline zrovector_small<l>& set(const CPPL_INT&, const comple&);
  template<CPPL_INT _l> inline friend std::ostream& operator<<(std::ostream&, const zrovector_small<_l>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// calc ////////
#ifndef _MSC_VER
  template<CPPL_INT _l> inline friend zcovector_small<_l> t(const zrovector_small<_l>&);
  template<CPPL_INT _l> inline friend comple nrm2(const zrovector_small<_l>&);
  template<CPPL_INT _l> inline friend CPPL_INT idamax(const zrovector_small<_l>&);
  template<CPPL_INT _l> inline friend comple damax(const zrovector_small<_l>&);
#endif//_MSC_VER
  
  //////// misc ////////
  inline zrovector_small<l>& zero();
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT L> inline zrovector_small<L>& operator= (const  zrovector_small<L>&);
  //////// += ////////
  template<CPPL_INT L> inline friend zrovector_small<L>& operator+=(zrovector_small<L>&, const zrovector_small<L>&);
  //////// -= ////////
  template<CPPL_INT L> inline friend zrovector_small<L>& operator-=(zrovector_small<L>&, const zrovector_small<L>&);
  //////// *= ////////
  template<CPPL_INT L> inline friend zrovector_small<L>& operator*=(zrovector_small<L>&, const             double&);
  template<CPPL_INT L> inline friend zrovector_small<L>& operator*=(zrovector_small<L>&, const             comple&);
  //////// /= ////////
  template<CPPL_INT L> inline friend zrovector_small<L>& operator/=(zrovector_small<L>&, const             double&);
  template<CPPL_INT L> inline friend zrovector_small<L>& operator/=(zrovector_small<L>&, const             comple&);
  //////// unary ////////
  template<CPPL_INT L> inline friend const zrovector_small<L>& operator+(const zrovector_small<L>&);
  template<CPPL_INT L> inline friend zrovector_small<L> operator-(const zrovector_small<L>&);
  //////// + ////////
  template<CPPL_INT L> inline friend zrovector_small<L> operator+(const zrovector_small<L>&, const zrovector_small<L>&);
  //////// - ////////
  template<CPPL_INT L> inline friend zrovector_small<L> operator-(const zrovector_small<L>&, const zrovector_small<L>&);
  //////// * ////////
  template<CPPL_INT L> inline friend             comple operator*(const zrovector_small<L>&, const zcovector_small<L>&);
  template<CPPL_INT M, CPPL_INT N> inline friend zrovector_small<N> operator*(const zrovector_small<M>&, const zgematrix_small<M,N>&);
  template<CPPL_INT L> inline friend zrovector_small<L> operator*(const zrovector_small<L>&, const dsymatrix_small<L>&);
  template<CPPL_INT L> inline friend zrovector_small<L> operator*(const zrovector_small<L>&, const             comple&);
  template<CPPL_INT L> inline friend zrovector_small<L> operator*(const zrovector_small<L>&, const             comple&);
  //////// / ////////
  template<CPPL_INT L> inline friend zrovector_small<L> operator/(const zrovector_small<L>&, const             double&);
  template<CPPL_INT L> inline friend zrovector_small<L> operator/(const zrovector_small<L>&, const             comple&);
  //////// comple ////////
  template<CPPL_INT L> inline friend zrovector_small<L> operator*(const             double&, const zrovector_small<L>&);
  template<CPPL_INT L> inline friend zrovector_small<L> operator*(const             comple&, const zrovector_small<L>&);
  //////// hadamerd ////////
  template<CPPL_INT L> inline friend zrovector_small<L>  hadamerd(const zrovector_small<L>&, const zrovector_small<L>&);  
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline comple  operator/(const zrovec2&, const zrovec2&);
inline zrovec3 operator/(const zrovec3&, const zrovec3&);
