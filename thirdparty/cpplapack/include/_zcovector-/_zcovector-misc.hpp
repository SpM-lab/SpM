//=============================================================================
/*! nullify all the vector data */
inline void _zcovector::nullify() const
{CPPL_VERBOSE_REPORT;
  l=0;
  array=NULL;
}

//=============================================================================
/*! destroy all the vector data */
inline void _zcovector::destroy() const
{CPPL_VERBOSE_REPORT;
  delete [] array;
  array=NULL;
}
