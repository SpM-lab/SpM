//=============================================================================
//! Samll Complex Double-precision Column Vector Class
template<CPPL_INT l> class zcovector_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  comple array[l];
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline zcovector_small();
  inline explicit zcovector_small(const zcovector&);
  inline zcovector_small(const comple&, const comple&);
  inline zcovector_small(const comple&, const comple&, const comple&);
  inline ~zcovector_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zcovector to_zcovector() const;
  
  //////// io ////////
  inline comple& operator()(const CPPL_INT&);
  inline comple  operator()(const CPPL_INT&) const;
  inline zcovector_small<l>& set(const CPPL_INT&, const comple&);
  template<CPPL_INT _l> inline friend std::ostream& operator<<(std::ostream&, const zcovector_small<_l>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// calc ////////
#ifndef _MSC_VER
  template<CPPL_INT _l> inline friend zrovector_small<_l> t(const zcovector_small<_l>&);
  template<CPPL_INT _l> inline friend comple nrm2(const zcovector_small<_l>&);
  template<CPPL_INT _l> inline friend CPPL_INT idamax(const zcovector_small<_l>&);
  template<CPPL_INT _l> inline friend comple damax(const zcovector_small<_l>&);
#endif//_MSC_VER
  
  //////// misc ////////
  inline zcovector_small<l>& zero();
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT L> inline zcovector_small<L>& operator= (const  zcovector_small<L>&);
  //////// += ////////
  template<CPPL_INT L> inline friend zcovector_small<L>& operator+=(zcovector_small<L>&, const zcovector_small<L>&);
  //////// -= ////////
  template<CPPL_INT L> inline friend zcovector_small<L>& operator-=(zcovector_small<L>&, const zcovector_small<L>&);
  //////// *= ////////
  template<CPPL_INT L> inline friend zcovector_small<L>& operator*=(zcovector_small<L>&, const             double&);
  template<CPPL_INT L> inline friend zcovector_small<L>& operator*=(zcovector_small<L>&, const             comple&);
  //////// /= ////////
  template<CPPL_INT L> inline friend zcovector_small<L>& operator/=(zcovector_small<L>&, const             double&);
  template<CPPL_INT L> inline friend zcovector_small<L>& operator/=(zcovector_small<L>&, const             comple&);
  //////// unary ////////
  template<CPPL_INT L> inline friend const zcovector_small<L>& operator+(const zcovector_small<L>&);
  template<CPPL_INT L> inline friend zcovector_small<L> operator-(const zcovector_small<L>&);
  //////// + ////////
  template<CPPL_INT L> inline friend zcovector_small<L> operator+(const zcovector_small<L>&, const zcovector_small<L>&);
  //////// - ////////
  template<CPPL_INT L> inline friend zcovector_small<L> operator-(const zcovector_small<L>&, const zcovector_small<L>&);
  //////// * ////////
  template<CPPL_INT M, CPPL_INT N> inline friend zgematrix_small<M,N> operator*(const zcovector_small<M>&, const zrovector_small<N>&);
  template<CPPL_INT L> inline friend zcovector_small<L> operator*(const zcovector_small<L>&, const             double&);
  template<CPPL_INT L> inline friend zcovector_small<L> operator*(const zcovector_small<L>&, const             comple&);
  //////// / ////////
  template<CPPL_INT L> inline friend zcovector_small<L> operator/(const zcovector_small<L>&, const             double&);
  template<CPPL_INT L> inline friend zcovector_small<L> operator/(const zcovector_small<L>&, const             comple&);
  //////// comple ////////
  template<CPPL_INT L> inline friend zcovector_small<L> operator*(const             double&, const zcovector_small<L>&);
  template<CPPL_INT L> inline friend zcovector_small<L> operator*(const             comple&, const zcovector_small<L>&);
  //////// hadamerd ////////
  template<CPPL_INT L> inline friend zcovector_small<L>  hadamerd(const zcovector_small<L>&, const zcovector_small<L>&);
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
//////// zcovec2 ////////
inline comple  operator/(const zcovec2&, const zcovec2&);

//////// zcovec3 ////////
inline zcovec3 operator/(const zcovec3&, const zcovec3&);
inline zcovec3 operator/=(zcovec3&, const zcovec3&);
