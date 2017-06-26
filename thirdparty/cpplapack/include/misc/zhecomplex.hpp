//=============================================================================
//! (DO NOT USE) Complex-double Class for Hermitian matrices
class zhecomplex : public comple
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT i, j;
  comple& v;
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline zhecomplex(const CPPL_INT&, const CPPL_INT&, comple&);
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  inline zhecomplex& operator=(const comple&);
  inline zhecomplex& operator+=(const comple&);
  inline zhecomplex& operator-=(const comple&);
  inline zhecomplex& operator*=(const comple&);
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! constructor */
inline zhecomplex::zhecomplex(const CPPL_INT& _i, const CPPL_INT& _j, comple& _v)
  : comple( _i < _j ? std::conj( _v ) : _v ), 
    v( _v )
{CPPL_VERBOSE_REPORT;
  i = _i;
  j = _j;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! operator= */
inline zhecomplex& zhecomplex::operator=(const comple& _v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i==j && std::fabs(_v.imag()) > DBL_MIN ){
    WARNING_REPORT;
    std::cerr << "Diagonal components of a hermitian matrix have to be real numbers." << std::endl
              << "Your input to the (" << i << "," << j << ") element was a complex number, " << _v << "." << std::endl;
  }
#endif//CPPL_DEBUG
  
  //comple::operator=( _v );
  //v = ( i < j ? std::conj( _v ) : _v );
  if(i>=j){
    v =_v;
  }
  else{//i<j
    v =std::conj(_v);
  }
  return *this;
}

//=============================================================================
/*! operator+= */
inline zhecomplex& zhecomplex::operator+=(const comple& _v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i==j && std::fabs(_v.imag()) > DBL_MIN ){
    WARNING_REPORT;
    std::cerr << "Diagonal components of a hermitian matrix have to be real numbers." << std::endl
              << "Your input to the (" << i << "," << j << ") element was a complex number, " << _v << "." << std::endl;
  }
#endif//CPPL_DEBUG
  
  if(i>=j){
    v +=_v;
  }
  else{//i<j
    v +=std::conj(_v);
  }
  return *this;
}

//=============================================================================
/*! operator-= */
inline zhecomplex& zhecomplex::operator-=(const comple& _v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i==j && std::fabs(_v.imag()) > DBL_MIN ){
    WARNING_REPORT;
    std::cerr << "Diagonal components of a hermitian matrix have to be real numbers." << std::endl
              << "Your input to the (" << i << "," << j << ") element was a complex number, " << _v << "." << std::endl;
  }
#endif//CPPL_DEBUG
  
  if(i>=j){
    v -=_v;
  }
  else{//i<j
    v -=std::conj(_v);
  }
  return *this;
}

//=============================================================================
/*! operator*= */
inline zhecomplex& zhecomplex::operator*=(const comple& _v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if( i==j && std::fabs(_v.imag()) > DBL_MIN ){
    WARNING_REPORT;
    std::cerr << "Diagonal components of a hermitian matrix have to be real numbers." << std::endl
              << "Your input to the (" << i << "," << j << ") element was a complex number, " << _v << "." << std::endl;
  }
#endif//CPPL_DEBUG
  
  if(i>=j){
    v *=_v;
  }
  else{//i<j
    v *=std::conj(_v);
  }
  return *this;
}
