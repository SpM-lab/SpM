//=============================================================================
/*! nullify all the matrix data */
inline void _zgsmatrix::nullify() const
{CPPL_VERBOSE_REPORT;
  m=0;
  n=0;
  data.clear();
  rows.clear();
  cols.clear();
}

//=============================================================================
/*! destroy all the matrix data */
inline void _zgsmatrix::destroy() const
{CPPL_VERBOSE_REPORT;
  data.clear();
  rows.clear();
  cols.clear();
}
