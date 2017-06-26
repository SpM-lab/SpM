//=============================================================================
/*! zhsmatrix*complex operator */
inline _zgsmatrix operator*(const zhsmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat( mat.to_zgsmatrix() );
  return newmat*d;
}

//=============================================================================
/*! zhsmatrix/complex operator */
inline _zgsmatrix operator/(const zhsmatrix& mat, const comple& d)
{CPPL_VERBOSE_REPORT;
  zgsmatrix newmat( mat.to_zgsmatrix() );
  return newmat/d;
}
