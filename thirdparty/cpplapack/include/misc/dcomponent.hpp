//=============================================================================
//! Component Class for Real Double-precision Sparse Matrix Classes
class dcomponent
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT i; //!< i index of the component
  CPPL_INT j; //!< j index of the component
  double v; //!< value of the component
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dcomponent(){ ; }
  inline dcomponent(const CPPL_INT& _i, const CPPL_INT& _j, const double& _v) :i(_i), j(_j), v(_v){ ; }
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  inline friend std::ostream& operator<<(std::ostream&, const dcomponent&);
  inline static bool ilt(const dcomponent&, const dcomponent&);
  inline static bool jlt(const dcomponent&, const dcomponent&);
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
inline std::ostream& operator<<(std::ostream& s, const dcomponent& c)
{CPPL_VERBOSE_REPORT;
  s << "(" << c.i << ", " << c.j << ",  " << c.v << ")" << std::flush;
  return s;
}

//=============================================================================
/*! lessthan function for i of dcomponent */
inline bool dcomponent::ilt(const dcomponent& a, const dcomponent& b)
{CPPL_VERBOSE_REPORT;
  return a.i < b.i;
}

//=============================================================================
/*! lessthan function for j of dcomponent */
inline bool dcomponent::jlt(const dcomponent& a, const dcomponent& b)
{CPPL_VERBOSE_REPORT;
  return a.j < b.j;
}
