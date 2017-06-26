//=============================================================================
//! Samll Real Double-precision Row Vector Class
template<CPPL_INT l> class drovector_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  double array[l]; //!< 1D array to store vector data
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline drovector_small();
  inline explicit drovector_small(const drovector&);
  inline drovector_small(const double&, const double&);
  inline drovector_small(const double&, const double&, const double&);
  inline drovector_small(const double&, const double&, const double&, const double&);
  inline ~drovector_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _drovector to_drovector() const;
  
  //////// io ////////
  inline double& operator()(const CPPL_INT&);
  inline double  operator()(const CPPL_INT&) const;
  inline drovector_small<l>& set(const CPPL_INT&, const double&);
  template<CPPL_INT _l> inline friend std::ostream& operator<<(std::ostream&, const drovector_small<_l>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// calc ////////
#ifndef _MSC_VER
  template<CPPL_INT _l> inline friend dcovector_small<_l> t(const drovector_small<_l>&);
  template<CPPL_INT _l> inline friend double nrm2(const drovector_small<_l>&);
  template<CPPL_INT _l> inline friend CPPL_INT idamax(const drovector_small<_l>&);
  template<CPPL_INT _l> inline friend double damax(const drovector_small<_l>&);
#endif//_MSC_VER
  
  //////// misc ////////
  inline drovector_small<l>& zero();
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT L> inline drovector_small<L>& operator= (const  drovector_small<L>&);
  //////// += ////////
  template<CPPL_INT L> inline friend drovector_small<L>& operator+=(drovector_small<L>&, const drovector_small<L>&);
  //////// -= ////////
  template<CPPL_INT L> inline friend drovector_small<L>& operator-=(drovector_small<L>&, const drovector_small<L>&);
  //////// *= ////////
  template<CPPL_INT L> inline friend drovector_small<L>& operator*=(drovector_small<L>&, const             double&);
  //////// /= ////////
  template<CPPL_INT L> inline friend drovector_small<L>& operator/=(drovector_small<L>&, const             double&);
  //////// unary ////////
  template<CPPL_INT L> inline friend const drovector_small<L>& operator+(const drovector_small<L>&);
  template<CPPL_INT L> inline friend drovector_small<L> operator-(const drovector_small<L>&);
  //////// + ////////
  template<CPPL_INT L> inline friend drovector_small<L> operator+(const drovector_small<L>&, const drovector_small<L>&);
  //////// - ////////
  template<CPPL_INT L> inline friend drovector_small<L> operator-(const drovector_small<L>&, const drovector_small<L>&);
  //////// * ////////
  template<CPPL_INT L> inline friend             double operator*(const drovector_small<L>&, const dcovector_small<L>&);
  template<CPPL_INT M, CPPL_INT N> inline friend drovector_small<N> operator*(const drovector_small<M>&, const dgematrix_small<M,N>&);
  template<CPPL_INT L> inline friend drovector_small<L> operator*(const drovector_small<L>&, const dsymatrix_small<L>&);
  template<CPPL_INT L> inline friend drovector_small<L> operator*(const drovector_small<L>&, const             double&);
  //////// / ////////
  template<CPPL_INT L> inline friend drovector_small<L> operator/(const drovector_small<L>&, const             double&);
  //////// double ////////
  template<CPPL_INT L> inline friend drovector_small<L> operator*(const             double&, const drovector_small<L>&);
  //////// hadamerd ////////
  template<CPPL_INT L> inline friend drovector_small<L>  hadamerd(const drovector_small<L>&, const drovector_small<L>&);
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

inline double  operator/(const drovec2&, const drovec2&);
inline drovec3 operator/(const drovec3&, const drovec3&);
