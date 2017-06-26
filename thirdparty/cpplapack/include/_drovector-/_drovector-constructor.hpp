//=============================================================================
/*! drovector constructor */
inline _drovector::_drovector()
{CPPL_VERBOSE_REPORT;
  l =0;
  cap =0;
  array =NULL;
}

//=============================================================================
/*! _drovector copy constructor */
inline _drovector::_drovector(const _drovector& vec)
{CPPL_VERBOSE_REPORT;
  l =vec.l;
  cap =vec.cap;
  array =vec.array;
  
  vec.nullify();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! _drovector destructor */
inline _drovector::~_drovector()
{CPPL_VERBOSE_REPORT;
  delete[] array; 
}
