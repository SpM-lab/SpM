//=============================================================================
//! Complex Double-precision Hermitian Sparse Matrix Class
class zhsmatrix
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT const& m; //!< matrix row size
  CPPL_INT n; //!< matrix column size
  std::vector<zcomponent> data; //!< matrix data
  std::vector< std::vector<CPPL_INT> > line; //!< vector of vector to store the entry information of component for each row and column
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline zhsmatrix();
  inline zhsmatrix(const zhsmatrix&);
  inline zhsmatrix(const _zhsmatrix&);
  inline zhsmatrix(const CPPL_INT&, const CPPL_INT=0);
  inline zhsmatrix(const char*);
  inline zhsmatrix(const zgematrix&);
  inline ~zhsmatrix(); //destructor
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zgematrix to_zgematrix() const;
  inline _zhematrix to_zhematrix() const;
  inline _zgsmatrix to_zgsmatrix() const;
  
  //////// io ////////
  inline comple operator()(const CPPL_INT&, const CPPL_INT&) const;
  inline zhecomplex operator()(const CPPL_INT&, const CPPL_INT&);
  inline zhsmatrix& put(const CPPL_INT&, const CPPL_INT&, const comple&);
  inline zhsmatrix& del(const CPPL_INT, const CPPL_INT); //<-- NOT (const CPPL_INT&)
  inline zhsmatrix& del(const CPPL_INT); //<-- NOT (const CPPL_INT&)
  inline friend std::ostream& operator<<(std::ostream&, const zhsmatrix&);
  inline void write(const char*) const;
  inline void read(const char*);

  //////// misc ////////
  inline void clear();
  inline zhsmatrix& zero();
  inline zhsmatrix& identity();
  inline void chsign();
  inline void copy(const zhsmatrix&);
  inline void shallow_copy(const _zhsmatrix&);
  inline zhsmatrix& resize(const CPPL_INT&, const CPPL_INT=0, const CPPL_INT=0);
  inline void stretch(const CPPL_INT&);
  inline void expand(const CPPL_INT&);
  inline bool isListed(const CPPL_INT&, const CPPL_INT&) const;
  inline CPPL_INT number(const CPPL_INT&, const CPPL_INT&) const;
  //inline void unique();
  inline void checkup();
  inline _zrovector row(const CPPL_INT&) const;
  inline _zcovector col(const CPPL_INT&) const;
  inline void diet(const double=DBL_MIN);
  inline friend void swap(zhsmatrix&, zhsmatrix&);
  inline friend _zhsmatrix _(zhsmatrix&);
  
  //////// calc ////////
  inline friend _zhsmatrix t(const zhsmatrix&);
  inline friend void idamax(CPPL_INT&, CPPL_INT&, const zhsmatrix&);
  inline friend comple damax(const zhsmatrix&);
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
  //////// = ////////
  inline zhsmatrix& operator=(const  zhsmatrix&);
  inline zhsmatrix& operator=(const _zhsmatrix&);
  
  //////// += ////////
  inline zhsmatrix& operator+=(const  zhsmatrix&);
  inline zhsmatrix& operator+=(const _zhsmatrix&);
  
  //////// -= ////////
  inline zhsmatrix& operator-=(const  zhsmatrix&);
  inline zhsmatrix& operator-=(const _zhsmatrix&);
  
  //////// *= ////////
  inline zhsmatrix& operator*=(const     double&);
  
  //////// /= ////////
  inline zhsmatrix& operator/=(const     double&);
  
  //////// unary ////////
  inline friend const zhsmatrix& operator+(const zhsmatrix&);
  inline friend _zhsmatrix operator-(const  zhsmatrix&);
  
  //////// + ////////
  inline friend _zgematrix operator+(const  zhsmatrix&, const  zgematrix&);
  inline friend _zgematrix operator+(const  zhsmatrix&, const _zgematrix&);
  inline friend _zgematrix operator+(const  zhsmatrix&, const  zhematrix&);
  inline friend _zgematrix operator+(const  zhsmatrix&, const _zhematrix&);
  inline friend _zgematrix operator+(const  zhsmatrix&, const  zgbmatrix&);
  inline friend _zgematrix operator+(const  zhsmatrix&, const _zgbmatrix&);
  inline friend _zgsmatrix operator+(const  zhsmatrix&, const  zgsmatrix&);
  inline friend _zgsmatrix operator+(const  zhsmatrix&, const _zgsmatrix&);
  inline friend _zhsmatrix operator+(const  zhsmatrix&, const  zhsmatrix&);
  inline friend _zhsmatrix operator+(const  zhsmatrix&, const _zhsmatrix&);
  
  //////// - ////////
  inline friend _zgematrix operator-(const  zhsmatrix&, const  zgematrix&);
  inline friend _zgematrix operator-(const  zhsmatrix&, const _zgematrix&);
  inline friend _zgematrix operator-(const  zhsmatrix&, const  zhematrix&);
  inline friend _zgematrix operator-(const  zhsmatrix&, const _zhematrix&);
  inline friend _zgematrix operator-(const  zhsmatrix&, const  zgbmatrix&);
  inline friend _zgematrix operator-(const  zhsmatrix&, const _zgbmatrix&);
  inline friend _zgsmatrix operator-(const  zhsmatrix&, const  zgsmatrix&);
  inline friend _zgsmatrix operator-(const  zhsmatrix&, const _zgsmatrix&);
  inline friend _zhsmatrix operator-(const  zhsmatrix&, const  zhsmatrix&);
  inline friend _zhsmatrix operator-(const  zhsmatrix&, const _zhsmatrix&);
  
  //////// * ////////
  inline friend _zcovector operator*(const  zhsmatrix&, const  zcovector&);
  inline friend _zcovector operator*(const  zhsmatrix&, const _zcovector&);
  inline friend _zgematrix operator*(const  zhsmatrix&, const  zgematrix&);
  inline friend _zgematrix operator*(const  zhsmatrix&, const _zgematrix&);  
  inline friend _zgematrix operator*(const  zhsmatrix&, const  zhematrix&);
  inline friend _zgematrix operator*(const  zhsmatrix&, const _zhematrix&);
  inline friend _zgematrix operator*(const  zhsmatrix&, const  zgbmatrix&);
  inline friend _zgematrix operator*(const  zhsmatrix&, const _zgbmatrix&);  
  inline friend _zgsmatrix operator*(const  zhsmatrix&, const  zgsmatrix&);
  inline friend _zgsmatrix operator*(const  zhsmatrix&, const _zgsmatrix&);
  inline friend _zgsmatrix operator*(const  zhsmatrix&, const  zhsmatrix&);
  inline friend _zgsmatrix operator*(const  zhsmatrix&, const _zhsmatrix&);
  inline friend _zhsmatrix operator*(const  zhsmatrix&, const     double&);
  inline friend _zgsmatrix operator*(const  zhsmatrix&, const     comple&);
  
  //////// / ////////
  inline friend _zhsmatrix operator/(const  zhsmatrix&, const     double&);
  inline friend _zgsmatrix operator/(const  zhsmatrix&, const     comple&);
  
  //////// double, comple ////////
  inline friend _zhsmatrix operator*(const     double&, const  zhsmatrix&);
  inline friend _zgsmatrix operator*(const     comple&, const  zhsmatrix&);
};
