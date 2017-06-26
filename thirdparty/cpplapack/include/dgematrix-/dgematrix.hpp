//=============================================================================
//! Real Double-precision General Dence Matrix Class
class dgematrix
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT m; //!< matrix row size
  CPPL_INT n; //!< matrix column size
  double* array; //!< 1D array to store matrix data
  double** darray; //!< array of pointers of column head addresses
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dgematrix();
  inline dgematrix(const dgematrix&);
  inline dgematrix(const _dgematrix&);
  inline dgematrix(const CPPL_INT&, const CPPL_INT&);
  inline dgematrix(const char*);
  inline ~dgematrix(); //destructor
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zgematrix to_zgematrix() const;
  
  //////// io ////////
  inline double& operator()(const CPPL_INT&, const CPPL_INT&);
  inline double operator()(const CPPL_INT&, const CPPL_INT&) const;
  inline dgematrix& set(const CPPL_INT&, const CPPL_INT&, const double&); //const;
  inline friend std::ostream& operator<<(std::ostream&, const dgematrix&);
  inline void write(const char*) const;
  inline void read(const char*);
  
  //////// misc ////////
  inline void clear();
  inline dgematrix& zero();
  inline dgematrix& identity();
  inline void chsign();
  inline void copy(const dgematrix&);
  inline void shallow_copy(const _dgematrix&);
  inline dgematrix& resize(const CPPL_INT&, const CPPL_INT&);
  inline _drovector row(const CPPL_INT&) const;
  inline _dcovector col(const CPPL_INT&) const;
  inline friend void swap(dgematrix&, dgematrix&);
  inline friend _dgematrix _(dgematrix&);
  
  //////// calc ////////
  inline friend _dgematrix t(const dgematrix&);
  inline friend _dgematrix i(const dgematrix&);
  inline friend void idamax(CPPL_INT&, CPPL_INT&, const dgematrix&);
  inline friend double damax(const dgematrix&);
  
  //////// lapack ////////
  inline CPPL_INT dgesv(dgematrix&);
  inline CPPL_INT dgesv(dcovector&);
  inline CPPL_INT dgels(dgematrix&);
  inline CPPL_INT dgels(dcovector&);
  inline CPPL_INT dgels(dgematrix&, drovector&);
  inline CPPL_INT dgels(dcovector&, double&);
  //inline CPPL_INT dgelss(dcovector&);
  inline CPPL_INT dgelss(dcovector&, dcovector&, CPPL_INT&, const double);
  inline CPPL_INT dgelss(dgematrix&, dcovector&, CPPL_INT&, const double);
  inline CPPL_INT dgelsd(dcovector&, dcovector&, CPPL_INT&, const double);
  //inline CPPL_INT dgelsd(dgematrix&, dcovector&, CPPL_INT&, const double);
  inline CPPL_INT dgeev(std::vector<double>&, std::vector<double>&);
  inline CPPL_INT dgeev(zcovector&);
  inline CPPL_INT dgeev(std::vector<double>&, std::vector<double>&, std::vector<dcovector>&, std::vector<dcovector>&);
  inline CPPL_INT dgeev(std::vector<double>&, std::vector<double>&, std::vector<drovector>&, std::vector<drovector>&);
  inline CPPL_INT dggev(dgematrix&, std::vector<double>&, std::vector<double>&);
  inline CPPL_INT dggev(dgematrix&, std::vector<double>&, std::vector<double>&, std::vector<dcovector>&, std::vector<dcovector>&);
  inline CPPL_INT dggev(dgematrix&, std::vector<double>&, std::vector<double>&, std::vector<drovector>&, std::vector<drovector>&);
  inline CPPL_INT dgesvd(dgbmatrix&);
  inline CPPL_INT dgesvd(dcovector&, dgematrix&, dgematrix&);
  inline CPPL_INT dgglse(dgematrix&, dcovector&, dcovector&, dcovector&);
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
  //////// = ////////
  inline dgematrix& operator=(const  dgematrix&);
  inline dgematrix& operator=(const _dgematrix&);
  
  //////// += ////////
  inline dgematrix& operator+=(const  dgematrix&);
  inline dgematrix& operator+=(const _dgematrix&);
  inline dgematrix& operator+=(const  dsymatrix&);
  inline dgematrix& operator+=(const _dsymatrix&);
  inline dgematrix& operator+=(const  dgbmatrix&);
  inline dgematrix& operator+=(const _dgbmatrix&);
  inline dgematrix& operator+=(const  dgsmatrix&);
  inline dgematrix& operator+=(const _dgsmatrix&);
  inline dgematrix& operator+=(const  dssmatrix&);
  inline dgematrix& operator+=(const _dssmatrix&);
  
  //////// -= ////////
  inline dgematrix& operator-=(const  dgematrix&);
  inline dgematrix& operator-=(const _dgematrix&);
  inline dgematrix& operator-=(const  dsymatrix&);
  inline dgematrix& operator-=(const _dsymatrix&);
  inline dgematrix& operator-=(const  dgbmatrix&);
  inline dgematrix& operator-=(const _dgbmatrix&);
  inline dgematrix& operator-=(const  dgsmatrix&);
  inline dgematrix& operator-=(const _dgsmatrix&);
  inline dgematrix& operator-=(const  dssmatrix&);
  inline dgematrix& operator-=(const _dssmatrix&);
  
  //////// *= ////////
  inline dgematrix& operator*=(const  dgematrix&);
  inline dgematrix& operator*=(const _dgematrix&);
  inline dgematrix& operator*=(const  dsymatrix&);
  inline dgematrix& operator*=(const _dsymatrix&);
  inline dgematrix& operator*=(const  dgbmatrix&);
  inline dgematrix& operator*=(const _dgbmatrix&);
  inline dgematrix& operator*=(const  dgsmatrix&);
  inline dgematrix& operator*=(const _dgsmatrix&);
  inline dgematrix& operator*=(const  dssmatrix&);
  inline dgematrix& operator*=(const _dssmatrix&);
  inline dgematrix& operator*=(const     double&);
  
  //////// /= ////////
  inline dgematrix& operator/=(const     double&);

  //////// unary ////////
  inline friend const dgematrix& operator+(const dgematrix&);
  inline friend _dgematrix operator-(const  dgematrix&);
  
  //////// + ////////
  inline friend _dgematrix operator+(const  dgematrix&, const  dgematrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const _dgematrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const  dsymatrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const _dsymatrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const  dgbmatrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const _dgbmatrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const  dgsmatrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const _dgsmatrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const  dssmatrix&);
  inline friend _dgematrix operator+(const  dgematrix&, const _dssmatrix&);
  
  //////// - ////////
  inline friend _dgematrix operator-(const  dgematrix&, const  dgematrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const _dgematrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const  dsymatrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const _dsymatrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const  dgbmatrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const _dgbmatrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const  dgsmatrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const _dgsmatrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const  dssmatrix&);
  inline friend _dgematrix operator-(const  dgematrix&, const _dssmatrix&);
  
  //////// * ////////
  inline friend _dcovector operator*(const  dgematrix&, const  dcovector&);
  inline friend _dcovector operator*(const  dgematrix&, const _dcovector&);
  inline friend _dgematrix operator*(const  dgematrix&, const  dgematrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const _dgematrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const  dsymatrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const _dsymatrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const  dgbmatrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const _dgbmatrix&);  
  inline friend _dgematrix operator*(const  dgematrix&, const  dgsmatrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const _dgsmatrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const  dssmatrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const _dssmatrix&);
  inline friend _dgematrix operator*(const  dgematrix&, const     double&);
  
  //////// / ////////
  inline friend _dgematrix operator/(const  dgematrix&, const     double&);
  
  //////// % ////////
  inline friend _drovector operator%(const  dgematrix&, const  dgematrix&);
  
  //////// double ////////
  inline friend _dgematrix operator*(const     double&, const  dgematrix&);
  
  //////// hadamard ////////
  inline friend _dgematrix hadamard(const dgematrix&, const dgematrix&);
};
