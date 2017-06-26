//=============================================================================
//! Real Double-precision General Sparse Matrix Class
class dgsmatrix
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT m; //!< matrix row size
  CPPL_INT n; //!< matrix column size
  std::vector<dcomponent> data; //!< matrix data
  std::vector< std::vector<CPPL_INT> > rows; //!< array of vector to store the entry information of component for each row
  std::vector< std::vector<CPPL_INT> > cols; //!< array of vector to store the entry information of component for each column
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dgsmatrix();
  inline dgsmatrix(const dgsmatrix&);
  inline dgsmatrix(const _dgsmatrix&);
  inline dgsmatrix(const CPPL_INT&, const CPPL_INT&, const CPPL_INT=0);
  inline dgsmatrix(const char*);
  inline ~dgsmatrix(); //destructor
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zgsmatrix to_zgsmatrix() const;
  inline _dgematrix to_dgematrix() const;
  inline  dgrmatrix to_dgrmatrix() const;
  
  //////// io ////////
  inline double operator()(const CPPL_INT&, const CPPL_INT&) const;
  inline double& operator()(const CPPL_INT&, const CPPL_INT&);
  inline dgsmatrix& put(const CPPL_INT&, const CPPL_INT&, const double&);
  inline dgsmatrix& del(const CPPL_INT, const CPPL_INT); //<-- NOT (const CPPL_INT&)
  inline dgsmatrix& del(const CPPL_INT); //<-- NOT (const CPPL_INT&)
  inline friend std::ostream& operator<<(std::ostream&, const dgsmatrix&);
  inline void write(const char*) const;
  inline void read(const char*);

  //////// misc ////////
  inline void clear();
  inline dgsmatrix& zero();
  inline dgsmatrix& identity();
  inline void chsign();
  inline void copy(const dgsmatrix&);
  inline void shallow_copy(const _dgsmatrix&);
  inline dgsmatrix& resize(const CPPL_INT&, const CPPL_INT&, const CPPL_INT=0, const CPPL_INT=0);
  inline void stretch(const CPPL_INT&, const CPPL_INT&);
  inline bool isListed(const CPPL_INT&, const CPPL_INT&) const;
  inline CPPL_INT number(const CPPL_INT&, const CPPL_INT&);
  inline void diet(const double=DBL_MIN);
  inline void checkup();
  inline _drovector row(const CPPL_INT&) const;
  inline _dcovector col(const CPPL_INT&) const;
  inline friend void swap(dgsmatrix&, dgsmatrix&);
  inline friend _dgsmatrix _(dgsmatrix&);
  
  //////// calc ////////
  inline friend _dgsmatrix t(const dgsmatrix&);
  inline friend void idamax(CPPL_INT&, CPPL_INT&, const dgsmatrix&);
  inline friend double damax(const dgsmatrix&);
  
  //////// pardiso ////////
  inline CPPL_INT pardiso(dcovector&) const;
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
  //////// = ////////
  inline dgsmatrix& operator=(const  dgsmatrix&);
  inline dgsmatrix& operator=(const _dgsmatrix&);
  
  //////// += ////////
  inline dgsmatrix& operator+=(const  dgsmatrix&);
  inline dgsmatrix& operator+=(const _dgsmatrix&);
  
  //////// -= ////////
  inline dgsmatrix& operator-=(const  dgsmatrix&);
  inline dgsmatrix& operator-=(const _dgsmatrix&);
  
  //////// *= ////////
  inline dgsmatrix& operator*=(const  dgsmatrix&);
  inline dgsmatrix& operator*=(const _dgsmatrix&);
  inline dgsmatrix& operator*=(const     double&);
  
  //////// /= ////////
  inline dgsmatrix& operator/=(const     double&);

  //////// unary ////////
  inline friend const dgsmatrix& operator+(const dgsmatrix&);
  inline friend _dgsmatrix operator-(const  dgsmatrix&);
  
  //////// + ////////
  inline friend _dgematrix operator+(const  dgsmatrix&, const  dgematrix&);
  inline friend _dgematrix operator+(const  dgsmatrix&, const _dgematrix&);
  inline friend _dgematrix operator+(const  dgsmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator+(const  dgsmatrix&, const _dsymatrix&);
  inline friend _dgematrix operator+(const  dgsmatrix&, const  dgbmatrix&);
  inline friend _dgematrix operator+(const  dgsmatrix&, const _dgbmatrix&);
  inline friend _dgsmatrix operator+(const  dgsmatrix&, const  dgsmatrix&);
  inline friend _dgsmatrix operator+(const  dgsmatrix&, const _dgsmatrix&);
  inline friend _dgsmatrix operator+(const  dgsmatrix&, const  dssmatrix&);
  inline friend _dgsmatrix operator+(const  dgsmatrix&, const _dssmatrix&);
  
  //////// - ////////
  inline friend _dgematrix operator-(const  dgsmatrix&, const  dgematrix&);
  inline friend _dgematrix operator-(const  dgsmatrix&, const _dgematrix&);
  inline friend _dgematrix operator-(const  dgsmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator-(const  dgsmatrix&, const _dsymatrix&);
  inline friend _dgematrix operator-(const  dgsmatrix&, const  dgbmatrix&);
  inline friend _dgematrix operator-(const  dgsmatrix&, const _dgbmatrix&);
  inline friend _dgsmatrix operator-(const  dgsmatrix&, const  dgsmatrix&);
  inline friend _dgsmatrix operator-(const  dgsmatrix&, const _dgsmatrix&);
  inline friend _dgsmatrix operator-(const  dgsmatrix&, const  dssmatrix&);
  inline friend _dgsmatrix operator-(const  dgsmatrix&, const _dssmatrix&);
  
  //////// * ////////
  inline friend _dcovector operator*(const  dgsmatrix&, const  dcovector&);
  inline friend _dcovector operator*(const  dgsmatrix&, const _dcovector&);
  inline friend _dgematrix operator*(const  dgsmatrix&, const  dgematrix&);
  inline friend _dgematrix operator*(const  dgsmatrix&, const _dgematrix&);  
  inline friend _dgematrix operator*(const  dgsmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator*(const  dgsmatrix&, const _dsymatrix&);
  inline friend _dgematrix operator*(const  dgsmatrix&, const  dgbmatrix&);
  inline friend _dgematrix operator*(const  dgsmatrix&, const _dgbmatrix&);  
  inline friend _dgsmatrix operator*(const  dgsmatrix&, const  dgsmatrix&);
  inline friend _dgsmatrix operator*(const  dgsmatrix&, const _dgsmatrix&);
  inline friend _dgsmatrix operator*(const  dgsmatrix&, const  dssmatrix&);
  inline friend _dgsmatrix operator*(const  dgsmatrix&, const _dssmatrix&);
  inline friend _dgsmatrix operator*(const  dgsmatrix&, const     double&);
  
  //////// / ////////
  inline friend _dgsmatrix operator/(const  dgsmatrix&, const     double&);

  //////// double ////////
  inline friend _dgsmatrix operator*(const     double&, const  dgsmatrix&);
};
