//=============================================================================
//! Real Double-precision General Band Matrix Class
class dgbmatrix
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT m; //!< matrix row size
  CPPL_INT n; //!< matrix column size
  CPPL_INT kl; //!< lower band width
  CPPL_INT ku; //!< upper band width
  double* array; //!< 1D array to store matrix data
  double** darray; //!< array of pointers of column head addresses
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dgbmatrix();
  inline dgbmatrix(const dgbmatrix&);
  inline dgbmatrix(const _dgbmatrix&);
  inline dgbmatrix(const CPPL_INT&, const CPPL_INT&, const CPPL_INT&, const CPPL_INT&);
  inline dgbmatrix(const char *);
  inline ~dgbmatrix(); //destructor
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zgbmatrix to_zgbmatrix() const;
  inline _dgematrix to_dgematrix() const;
  
  //////// io ////////
  inline double& operator()(const CPPL_INT&, const CPPL_INT&);
  inline double operator()(const CPPL_INT&, const CPPL_INT&) const;
  inline dgbmatrix& set(const CPPL_INT&, const CPPL_INT&, const double&); //const;
  inline friend std::ostream& operator<<(std::ostream&, const dgbmatrix&);
  inline void write(const char*) const;
  inline void read(const char*);
  
  //////// misc ////////
  inline void clear();
  inline dgbmatrix& zero();
  inline dgbmatrix& identity();
  inline void chsign();
  inline void copy(const dgbmatrix&);
  inline void shallow_copy(const _dgbmatrix&);
  inline dgbmatrix& resize(const CPPL_INT&, const CPPL_INT&, const CPPL_INT&, const CPPL_INT&);
  inline _drovector row(const CPPL_INT&) const;
  inline _dcovector col(const CPPL_INT&) const;
  inline friend void swap(dgbmatrix&, dgbmatrix&);
  inline friend _dgbmatrix _(dgbmatrix&);
  
  //////// calc ////////
  inline friend _dgbmatrix t(const dgbmatrix&);
  inline friend _dgematrix i(const dgbmatrix&);
  
  //////// lapack ////////
  inline CPPL_INT dgbsv(dgematrix&);
  inline CPPL_INT dgbsv(dcovector&);
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
  //////// = ////////
  inline dgbmatrix& operator=(const  dgbmatrix&);
  inline dgbmatrix& operator=(const _dgbmatrix&);
  
  //////// += ////////
  inline dgbmatrix& operator+=(const  dgbmatrix&);
  inline dgbmatrix& operator+=(const _dgbmatrix&);
  
  //////// -= ////////
  inline dgbmatrix& operator-=(const  dgbmatrix&);
  inline dgbmatrix& operator-=(const _dgbmatrix&);
  
  //////// *= ////////
  inline dgbmatrix& operator*=(const  dgbmatrix&);
  inline dgbmatrix& operator*=(const _dgbmatrix&);
  inline dgbmatrix& operator*=(const     double&);
  
  //////// /= ////////
  inline dgbmatrix& operator/=(const     double&);
  
  //////// unary ////////
  inline friend const dgbmatrix& operator+(const dgbmatrix&);
  inline friend _dgbmatrix operator-(const dgbmatrix&);
  
  //////// + ////////
  inline friend _dgematrix operator+(const  dgbmatrix&, const  dgematrix&);
  inline friend _dgematrix operator+(const  dgbmatrix&, const _dgematrix&);
  inline friend _dgematrix operator+(const  dgbmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator+(const  dgbmatrix&, const _dsymatrix&);
  inline friend _dgbmatrix operator+(const  dgbmatrix&, const  dgbmatrix&);
  inline friend _dgbmatrix operator+(const  dgbmatrix&, const _dgbmatrix&);
  inline friend _dgematrix operator+(const  dgbmatrix&, const  dgsmatrix&);
  inline friend _dgematrix operator+(const  dgbmatrix&, const _dgsmatrix&);
  inline friend _dgematrix operator+(const  dgbmatrix&, const  dssmatrix&);
  inline friend _dgematrix operator+(const  dgbmatrix&, const _dssmatrix&);
  
  //////// - ////////
  inline friend _dgematrix operator-(const  dgbmatrix&, const  dgematrix&);
  inline friend _dgematrix operator-(const  dgbmatrix&, const _dgematrix&);
  inline friend _dgematrix operator-(const  dgbmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator-(const  dgbmatrix&, const _dsymatrix&);
  inline friend _dgbmatrix operator-(const  dgbmatrix&, const  dgbmatrix&);
  inline friend _dgbmatrix operator-(const  dgbmatrix&, const _dgbmatrix&);
  inline friend _dgematrix operator-(const  dgbmatrix&, const  dgsmatrix&);
  inline friend _dgematrix operator-(const  dgbmatrix&, const _dgsmatrix&);
  inline friend _dgematrix operator-(const  dgbmatrix&, const  dssmatrix&);
  inline friend _dgematrix operator-(const  dgbmatrix&, const _dssmatrix&);
  
  //////// * ////////
  inline friend _dcovector operator*(const  dgbmatrix&, const  dcovector&);
  inline friend _dcovector operator*(const  dgbmatrix&, const _dcovector&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const  dgematrix&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const _dgematrix&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const  dsymatrix&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const _dsymatrix&);
  inline friend _dgbmatrix operator*(const  dgbmatrix&, const  dgbmatrix&);
  inline friend _dgbmatrix operator*(const  dgbmatrix&, const _dgbmatrix&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const  dgsmatrix&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const _dgsmatrix&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const  dssmatrix&);
  inline friend _dgematrix operator*(const  dgbmatrix&, const _dssmatrix&);
  inline friend _dgbmatrix operator*(const  dgbmatrix&, const     double&);
  
  //////// / ////////
  inline friend _dgbmatrix operator/(const  dgbmatrix&, const     double&);
  
  //////// double ////////
  inline friend _dgbmatrix operator*(const     double&, const  dgbmatrix&);
};
