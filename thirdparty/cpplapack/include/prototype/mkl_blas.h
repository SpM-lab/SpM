#define  drot_   DROT
#define  dswap_  DSWAP
#define  dscal_  DSCAL
#define  dcopy_  DCOPY
#define  daxpy_  DAXPY
#define  ddot_   DDOT
#define  dnrm2_  DNRM2
#define  dasum_  DASUM
#define  idamax_ IDAMAX
  
#define  dgemv_  DGEMV
#define  dgbmv_  DGBMV
#define  dsymv_  DSYMV
#define  dsbmv_  DSBMV
#define  dspmv_  DSPMV
#define  dtrmv_  DTRMV
#define  dtbmv_  DTBMV
#define  dtpmv_  DTPMV
#define  dtrsv_  DTRSV
#define  dtbsv_  DTBSV
#define  dtpsv_  DTPSV
#define  dger_   DGER
#define  dsyr_   DSYR
#define  dspr_   DSPR
#define  dsyr2_  DSYR2
#define  dspr2_  DSPR2

#define  dgemm_  DGEMM
#define  dsymm_  DSYMM
#define  dsyrk_  DSYRK
#define  dsyr2k_ DSYR2K
#define  dtrmm_  DTRMM
#define  dtrsm_  DTRSM

#define  zdrot_  ZDROT
#define  zswap_  ZSWAP
#define  zdscal_ ZDSCAL
#define  zscal_  ZSCAL
#define  zcopy_  ZCOPY
#define  zaxpy_  ZAXPY
#define  zdotu_  ZDOTU
#define  zdotc_  ZDOTC
#define  dznrm2_ DZNRM2
#define  dzasum_ DZASUM
#define  izamax_ IZAMAX
                  
#define  zgemv_  ZGEMV
#define  zgbmv_  ZGBMV
#define  zhemv_  ZHEMV
#define  zhbmv_  ZHBMV
#define  zhpmv_  ZHPMV
#define  ztrmv_  ZTRMV
#define  ztbmv_  ZTBMV
#define  ztpmv_  ZTPMV
#define  ztrsv_  ZTRSV
#define  ztbsv_  ZTBSV
#define  ztpsv_  ZTPSV
#define  zgeru_  ZGERU
#define  zgerc_  ZGERC
#define  zher_   ZHER
#define  zhpr_   ZHPR
#define  zher2_  ZHER2
#define  zhpr2_  ZHPR2
  
#define  zgemm_  ZGEMM
#define  zsymm_  ZSYMM
#define  zhemm_  ZHEMM
#define  zsyrk_  ZSYRK
#define  zherk_  ZHERK
#define  zsyr2k_ ZSYR2K
#define  zher2k_ ZHER2K
#define  ztrmm_  ZTRMM
#define  ztrsm_  ZTRSM

//////// mkl.h does not have legal zdotu_ and zdotc_ . Why??? ////////
inline MKL_Complex16 zdotu_( const CPPL_INT *N, MKL_Complex16 *x, const CPPL_INT *incx, const MKL_Complex16 *y, const CPPL_INT *incy )
{
  MKL_Complex16 v;
  zdotu_(&v, N, x, incx, y, incy);
  return v;
}

inline MKL_Complex16 zdotc_( const CPPL_INT *N, MKL_Complex16 *x, const CPPL_INT *incx, const MKL_Complex16 *y, const CPPL_INT *incy )
{
  MKL_Complex16 v;
  zdotc_(&v, N, x, incx, y, incy);
  return v;
}
