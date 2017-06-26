//=============================================================================
/*! _zrovector*zcovector operator */
inline comple operator*(const _zrovector& rovec, const zcovector& covec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(rovec.l!=covec.l){
    ERROR_REPORT;
    std::cerr << "These two vectors can not make a product." << std::endl
              << "Your input was (" << rovec.l << ") * (" << covec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  CPPL_INT inc =1;
  
  comple val =zdotu_( &rovec.l, rovec.array, &inc, covec.array, &inc );
  
  rovec.destroy();
  return val;
}
