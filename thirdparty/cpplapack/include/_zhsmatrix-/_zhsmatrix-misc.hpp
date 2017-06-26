//=============================================================================
/*! nullify all the matrix data */
inline void _zhsmatrix::nullify() const
{CPPL_VERBOSE_REPORT;
  n=0;
  data.clear();
  line.clear();
}

//=============================================================================
/*! destroy all the matrix data */
inline void _zhsmatrix::destroy() const
{CPPL_VERBOSE_REPORT;
  data.clear();
  line.clear();
}
