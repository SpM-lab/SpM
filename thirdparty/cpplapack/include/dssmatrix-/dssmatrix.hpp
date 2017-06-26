//=============================================================================
//! Real Double-precision Symmetric Sparse Matrix Class
class dssmatrix
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT const& m; //!< matrix row size
  CPPL_INT n; //!< matrix column size
  std::vector<dcomponent> data; //!< matrix data
  std::vector< std::vector<CPPL_INT> > line; //!< vector of vector to store the entry information of component for each row and column
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dssmatrix();
  inline dssmatrix(const dssmatrix&);
  inline dssmatrix(const _dssmatrix&);
  inline dssmatrix(const CPPL_INT&, const CPPL_INT=0);
  inline dssmatrix(const char*);
  inline ~dssmatrix(); //destructor
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zhsmatrix to_zhsmatrix() const;
  inline _dgematrix to_dgematrix() const;
  inline _dsymatrix to_dsymatrix() const;
  inline _dgsmatrix to_dgsmatrix() const;
  
  //////// io ////////
  inline double operator()(const CPPL_INT&, const CPPL_INT&) const;
  inline double& operator()(const CPPL_INT&, const CPPL_INT&);
  inline dssmatrix& put(const CPPL_INT&, const CPPL_INT&, const double&);
  inline dssmatrix& del(const CPPL_INT, const CPPL_INT); //<-- NOT (const CPPL_INT&)
  inline dssmatrix& del(const CPPL_INT); //<-- NOT (const CPPL_INT&)
  inline friend std::ostream& operator<<(std::ostream&, const dssmatrix&);
  inline void write(const char*) const;
  inline void read(const char*);

  //////// misc ////////
  inline void clear();
  inline dssmatrix& zero();
  inline dssmatrix& identity();
  inline void chsign();
  inline void copy(const dssmatrix&);
  inline void shallow_copy(const _dssmatrix&);
  inline dssmatrix& resize(const CPPL_INT&, const CPPL_INT=0, const CPPL_INT=0);
  inline void stretch(const CPPL_INT&);
  inline bool isListed(const CPPL_INT&, const CPPL_INT&) const;
  inline CPPL_INT number(const CPPL_INT&, const CPPL_INT&) const;
  inline _drovector row(const CPPL_INT&) const;
  inline _dcovector col(const CPPL_INT&) const;
  inline void diet(const double=DBL_MIN);
  inline CPPL_INT diag_front();
  inline void reorder(const bool=0);
  inline void rebuild();
  inline void checkup();
  inline friend void swap(dssmatrix&, dssmatrix&);
  inline friend _dssmatrix _(dssmatrix&);
  
  //////// calc ////////
  inline friend _dssmatrix t(const dssmatrix&);
  inline friend void idamax(CPPL_INT&, CPPL_INT&, const dssmatrix&);
  inline friend double damax(const dssmatrix&);
  
  //////// pardiso ////////
  inline CPPL_INT pardiso_definite(dcovector&) const;
  inline CPPL_INT pardiso_indefinite(dcovector&) const;
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
  //////// = ////////
  inline dssmatrix& operator=(const  dssmatrix&);
  inline dssmatrix& operator=(const _dssmatrix&);
  
  //////// += ////////
  inline dssmatrix& operator+=(const  dssmatrix&);
  inline dssmatrix& operator+=(const _dssmatrix&);
  
  //////// -= ////////
  inline dssmatrix& operator-=(const  dssmatrix&);
  inline dssmatrix& operator-=(const _dssmatrix&);
  
  //////// *= ////////
  inline dssmatrix& operator*=(const     double&);
  
  //////// /= ////////
  inline dssmatrix& operator/=(const     double&);
  
  //////// unary ////////
  inline friend const dssmatrix& operator+(const dssmatrix&);
  inline friend _dssmatrix operator-(const  dssmatrix&);
  
  //////// + ////////
  inline friend _dgematrix operator+(const  dssmatrix&, const  dgematrix&);
  inline friend _dgematrix operator+(const  dssmatrix&, const _dgematrix&);
  inline friend _dgematrix operator+(const  dssmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator+(const  dssmatrix&, const _dsymatrix&);
  inline friend _dgematrix operator+(const  dssmatrix&, const  dgbmatrix&);
  inline friend _dgematrix operator+(const  dssmatrix&, const _dgbmatrix&);
  inline friend _dgsmatrix operator+(const  dssmatrix&, const  dgsmatrix&);
  inline friend _dgsmatrix operator+(const  dssmatrix&, const _dgsmatrix&);
  inline friend _dssmatrix operator+(const  dssmatrix&, const  dssmatrix&);
  inline friend _dssmatrix operator+(const  dssmatrix&, const _dssmatrix&);
  
  //////// - ////////
  inline friend _dgematrix operator-(const  dssmatrix&, const  dgematrix&);
  inline friend _dgematrix operator-(const  dssmatrix&, const _dgematrix&);
  inline friend _dgematrix operator-(const  dssmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator-(const  dssmatrix&, const _dsymatrix&);
  inline friend _dgematrix operator-(const  dssmatrix&, const  dgbmatrix&);
  inline friend _dgematrix operator-(const  dssmatrix&, const _dgbmatrix&);
  inline friend _dgsmatrix operator-(const  dssmatrix&, const  dgsmatrix&);
  inline friend _dgsmatrix operator-(const  dssmatrix&, const _dgsmatrix&);
  inline friend _dssmatrix operator-(const  dssmatrix&, const  dssmatrix&);
  inline friend _dssmatrix operator-(const  dssmatrix&, const _dssmatrix&);
  
  //////// * ////////
  inline friend _dcovector operator*(const  dssmatrix&, const  dcovector&);
  inline friend _dcovector operator*(const  dssmatrix&, const _dcovector&);
  inline friend _dgematrix operator*(const  dssmatrix&, const  dgematrix&);
  inline friend _dgematrix operator*(const  dssmatrix&, const _dgematrix&);  
  inline friend _dgematrix operator*(const  dssmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator*(const  dssmatrix&, const _dsymatrix&);
  inline friend _dgematrix operator*(const  dssmatrix&, const  dgbmatrix&);
  inline friend _dgematrix operator*(const  dssmatrix&, const _dgbmatrix&);  
  inline friend _dgsmatrix operator*(const  dssmatrix&, const  dgsmatrix&);
  inline friend _dgsmatrix operator*(const  dssmatrix&, const _dgsmatrix&);  
  inline friend _dgsmatrix operator*(const  dssmatrix&, const  dssmatrix&);
  inline friend _dgsmatrix operator*(const  dssmatrix&, const _dssmatrix&);
  inline friend _dssmatrix operator*(const  dssmatrix&, const     double&);
  
  //////// / ////////
  inline friend _dssmatrix operator/(const  dssmatrix&, const     double&);
  
  //////// double ////////
  inline friend _dssmatrix operator*(const     double&, const  dssmatrix&);
};
