//=============================================================================
/*! convert to _zgematrix */
inline _zgematrix zgsmatrix::to_zgematrix() const
{CPPL_VERBOSE_REPORT;
  zgematrix newmat( zgematrix(m,n).zero() );
  
  const size_t data_size =data.size();
  for(size_t c=0; c<data_size; c++){
    const zcomponent& z =data[c];
    newmat(z.i,z.j) =z.v;
  }
  
  return _(newmat);
}
