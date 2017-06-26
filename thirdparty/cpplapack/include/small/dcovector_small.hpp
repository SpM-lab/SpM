//=============================================================================
//! Samll Real Double-precision Column Vector Class
template<CPPL_INT l> class dcovector_small
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  double array[l];
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dcovector_small();
  inline explicit dcovector_small(const dcovector&);
  inline dcovector_small(const double&, const double&);
  inline dcovector_small(const double&, const double&, const double&);
  inline dcovector_small(const double&, const double&, const double&, const double&);
  inline ~dcovector_small();
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _dcovector to_dcovector() const;
  
  //////// io ////////
  inline double& operator()(const CPPL_INT&);
  inline double  operator()(const CPPL_INT&) const;
  inline dcovector_small<l>& set(const CPPL_INT&, const double&);
  template<CPPL_INT _l> inline friend std::ostream& operator<<(std::ostream&, const dcovector_small<_l>&);
  inline void read(const char* filename);
  inline void write(const char* filename) const;
  
  //////// calc ////////
#ifndef _MSC_VER
  template<CPPL_INT _l> inline friend drovector_small<_l> t(const dcovector_small<_l>&);
  template<CPPL_INT _l> inline friend double nrm2(const dcovector_small<_l>&);
  template<CPPL_INT _l> inline friend CPPL_INT idamax(const dcovector_small<_l>&);
  template<CPPL_INT _l> inline friend double damax(const dcovector_small<_l>&);
#endif//_MSC_VER
  
  //////// misc ////////
  inline dcovector_small<l>& zero();
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
#ifndef _MSC_VER
  //////// = ////////
  template<CPPL_INT L> inline dcovector_small<L>& operator= (const  dcovector_small<L>&);
  //////// += ////////
  template<CPPL_INT L> inline friend dcovector_small<L>& operator+=(dcovector_small<L>&, const dcovector_small<L>&);
  //////// -= ////////
  template<CPPL_INT L> inline friend dcovector_small<L>& operator-=(dcovector_small<L>&, const dcovector_small<L>&);
  //////// *= ////////
  template<CPPL_INT L> inline friend dcovector_small<L>& operator*=(dcovector_small<L>&, const             double&);
  //////// /= ////////
  template<CPPL_INT L> inline friend dcovector_small<L>& operator/=(dcovector_small<L>&, const             double&);
  //////// unary ////////
  template<CPPL_INT L> inline friend const dcovector_small<L>& operator+(const dcovector_small<L>&);
  template<CPPL_INT L> inline friend dcovector_small<L> operator-(const dcovector_small<L>&);
  //////// + ////////
  template<CPPL_INT L> inline friend dcovector_small<L> operator+(const dcovector_small<L>&, const dcovector_small<L>&);
  //////// - ////////
  template<CPPL_INT L> inline friend dcovector_small<L> operator-(const dcovector_small<L>&, const dcovector_small<L>&);
  //////// * ////////
  template<CPPL_INT M, CPPL_INT N> inline friend dgematrix_small<M,N> operator*(const dcovector_small<M>&, const drovector_small<N>&);
  template<CPPL_INT L> inline friend dcovector_small<L> operator*(const dcovector_small<L>&, const             double&);
  //////// / ////////
  template<CPPL_INT L> inline friend dcovector_small<L> operator/(const dcovector_small<L>&, const             double&);
  //////// double ////////
  template<CPPL_INT L> inline friend dcovector_small<L> operator*(const             double&, const dcovector_small<L>&);
  //////// hadamerd ////////
  template<CPPL_INT L> inline friend dcovector_small<L>  hadamerd(const dcovector_small<L>&, const dcovector_small<L>&);
#endif//_MSC_VER
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
//////// dcovec2 ////////
inline double  operator/(const dcovec2&, const dcovec2&);
//inline dquater vr2q(const dcovec2&, const double&);
//inline dquater vt2q(const dcovec2&, const double&);
inline dcovec2 rotate(const dcovec2&, const double&);

//////// dcovec3 ////////
inline dcovec3 operator/(const dcovec3&, const dcovec3&);
inline dcovec3 operator/=(dcovec3&, const dcovec3&);
inline dquater vr2q(const dcovec3&, const double&);
inline dquater vt2q(const dcovec3&, const double&);
inline dcovec3 rotate(const dcovec3&, const dquater&);

//////// dquater ////////
inline dquater conj(const dquater&);
inline dcovec3 imag(const dquater&);
inline dquater inv(const dquater&);
inline dquater operator*(const dquater&, const dquater&);
inline dquater operator/(const dquater&, const dquater&);
inline dquater operator*=(dquater&, const dquater&);
inline dquater operator/=(dquater&, const dquater&);
inline dcovec3 q2vt(const dquater&);
inline dgemat3 q2m(const dquater&);
