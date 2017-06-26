//=============================================================================
//! Component Class for Complex Double-precision Sparse Matrix Classes
class zcomponent
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT i; //!< i index of the component
  CPPL_INT j; //!< j index of the component
  comple v; //!< value of the component
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline zcomponent(){ ; }
  inline zcomponent(const CPPL_INT& _i, const CPPL_INT& _j, const comple& _v) :i(_i), j(_j), v(_v){ ; }
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  inline friend std::ostream& operator<<(std::ostream&, const zcomponent&);
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const zcomponent& c)
{CPPL_VERBOSE_REPORT;
  s << "(" << c.i << ", " << c.j << ",  " << c.v << ")" << std::flush;
  return s;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! lessthan function for i of zcomponent */
inline bool ilt(const zcomponent& a, const zcomponent& b)
{CPPL_VERBOSE_REPORT;
  return a.i < b.i;
}

//=============================================================================
/*! lessthan function for j of zcomponent */
inline bool jlt(const zcomponent& a, const zcomponent& b)
{CPPL_VERBOSE_REPORT;
  return a.j < b.j;
}
