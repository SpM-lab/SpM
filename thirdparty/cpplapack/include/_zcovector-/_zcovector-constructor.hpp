//=============================================================================
/*! _zcovector constructor */
inline _zcovector::_zcovector()
{CPPL_VERBOSE_REPORT;
  l =0;
  array =NULL;
}

//=============================================================================
/*! _zcovector copy constructor */
inline _zcovector::_zcovector(const _zcovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  array =vec.array;
  
  vec.nullify();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! _zcovector destructor */
inline _zcovector::~_zcovector()
{CPPL_VERBOSE_REPORT;
  delete [] array;
}
