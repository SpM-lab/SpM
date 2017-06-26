//=============================================================================
/*! _dcovector constructor */
inline _dcovector::_dcovector()
{CPPL_VERBOSE_REPORT;CPPL_VERBOSE_REPORT;
  l =0;
  cap =0;
  array =NULL;
}

//=============================================================================
/*! _dcovector copy constructor */
inline _dcovector::_dcovector(const _dcovector& vec)
{CPPL_VERBOSE_REPORT;CPPL_VERBOSE_REPORT;
  l =vec.l;
  cap =vec.cap;
  array =vec.array;
  
  vec.nullify();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! _dcovector destructor */
inline _dcovector::~_dcovector()
{CPPL_VERBOSE_REPORT;CPPL_VERBOSE_REPORT;
  delete[] array;
}
