//=============================================================================
//! Real Double-precision General Compressed Sparse Row (CSR) Matrix Class
class dgrmatrix
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT m; //!< matrix row size
  CPPL_INT n; //!< matrix column size
  std::vector<double> a; //!< matrix component values
  std::vector<CPPL_INT> ia; //!< rowIndex (NOT zero-based BUT one-based indexing)
  std::vector<CPPL_INT> ja; //!< columns (NOT zero-based BUT one-based indexing)
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dgrmatrix();
  inline dgrmatrix(const CPPL_INT&, const CPPL_INT&);
  inline dgrmatrix(const dgrmatrix&);
  inline dgrmatrix(const char*);
  inline ~dgrmatrix(); //destructor
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  //inline _zgsmatrix to_zgsmatrix() const;
  inline _dgematrix to_dgematrix() const;
  
  //////// io ////////
  inline double operator()(const CPPL_INT&, const CPPL_INT&) const;
  inline double& operator()(const CPPL_INT&, const CPPL_INT&);
  inline friend std::ostream& operator<<(std::ostream&, const dgrmatrix&);
  inline void write(const char*) const;
  inline void read(const char*);
  
  //////// misc ////////
  inline void clear();
  inline dgrmatrix& zero();
  inline void copy(const dgrmatrix&);
  inline bool isListed(const CPPL_INT&, const CPPL_INT&) const;
  inline void checkup();
  //inline _drovector row(const CPPL_INT&) const;
  //inline _dcovector col(const CPPL_INT&) const;
  inline friend void swap(dgrmatrix&, dgrmatrix&);
  
  //////// calc ////////
  //inline friend void idamax(CPPL_INT&, CPPL_INT&, const dgrmatrix&);
  inline friend double damax(const dgrmatrix&);
  
  //////// MKL ////////
  inline CPPL_INT pardiso(dcovector&) const;
  inline CPPL_INT dfgmres(dcovector&, const double) const;
  inline CPPL_INT ilut_dfgmres(dcovector&, const int, const double) const;
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
  //////// = ////////
  inline dgrmatrix& operator=(const  dgrmatrix&);
  
  //////// ?= ////////
  inline dgrmatrix& operator*=(const     double&);
  inline dgrmatrix& operator/=(const     double&);  
  
  //////// * ////////
  inline friend  dgrmatrix operator*(const  dgrmatrix&, const     double&);
  inline friend  dgrmatrix operator*(const     double&, const  dgrmatrix&);
  inline friend _dcovector operator*(const  dgrmatrix&, const  dcovector&);
  inline friend _dcovector operator*(const  dgrmatrix&, const _dcovector&);
  
  //////// / ////////
  inline friend dgrmatrix operator/(const  dgrmatrix&, const     double&);
};
