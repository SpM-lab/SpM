extern "C" {
  // Level 1 BLAS
  /* Vector rotation: x := c*x[i] + s*y[i], y[i] := - s*x[i] + c*y[i] */
  void zdrot_( const CPPL_INT *N, std::complex<double> *x, const CPPL_INT *incx,
               std::complex<double> *y, const CPPL_INT *incy,
               const double *c, const double *s );
  /* x <--> y */
  void zswap_( const CPPL_INT *N, std::complex<double> *x, const CPPL_INT *incx,
               std::complex<double> *y, const CPPL_INT *incy );
  /* x := alpha * x  (alpha: double) */
  void zdscal_( const CPPL_INT *N, const double *alpha,
                std::complex<double> *x, const CPPL_INT *incx );
  /* x := alpha * x  (alpha: double complex)  */
  void zscal_( const CPPL_INT *N, const std::complex<double> *alpha,
               std::complex<double> *x, const CPPL_INT *incx );
  /* y := x */
  void zcopy_( const CPPL_INT *N, const std::complex<double> *x, const CPPL_INT *incx,
               std::complex<double> *y, const CPPL_INT *incy );
  /* y := alpha * x + y */
  void zaxpy_( const CPPL_INT *N, const std::complex<double> *alpha,
               const std::complex<double> *x, const CPPL_INT *incx,
               std::complex<double> *y, const CPPL_INT *incy );
  /* return x^T y */
  std::complex<double> zdotu_( const CPPL_INT *N, const std::complex<double> *x,
                               const CPPL_INT *incx, const std::complex<double> *y,
                               const CPPL_INT *incy );
  /* return conjg(x)^T y */
  std::complex<double> zdotc_( const CPPL_INT *N, const std::complex<double> *x,
                               const CPPL_INT *incx, const std::complex<double> *y,
                               const CPPL_INT *incy );
  /* return N2 norm */
  double dznrm2_( const CPPL_INT *N, const std::complex<double> *x,
                  const CPPL_INT *incx );
  /* return sum of abs(x[i]) */
  double dzasum_( const CPPL_INT *N, const std::complex<double> *x,
                  const CPPL_INT *incx );
  /* return the index of element having the largest absolute value */
  CPPL_INT izamax_( const CPPL_INT *N, const std::complex<double> *x, const CPPL_INT *incx );
  
  // Level 2 BLAS
  /* y := alpha * A * x + beta * y for General M by N Matrix */
  void zgemv_( const char *trans, const CPPL_INT *M, const CPPL_INT *N,
               const std::complex<double> *alpha, const std::complex<double> *a,
               const CPPL_INT *lda, const std::complex<double> *x, const CPPL_INT *incx,
               const std::complex<double> *beta, std::complex<double> *y,
               const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for General M by N Band Matrix */
  void zgbmv_( const char *trans, const CPPL_INT *M, const CPPL_INT *N,
               const CPPL_INT *KL, const CPPL_INT *KU,
               const std::complex<double> *alpha, const std::complex<double> *a,
               const CPPL_INT *lda, const std::complex<double> *x, const CPPL_INT *incx,
               const std::complex<double> *beta, std::complex<double> *y,
               const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for Hermitian Matrix */
  void zhemv_( const char *uplo, const CPPL_INT *N,
               const std::complex<double> *alpha, const std::complex<double> *a,
               const CPPL_INT *lda, const std::complex<double> *x, const CPPL_INT *incx,
               const std::complex<double> *beta, std::complex<double> *y,
               const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for Hermitian Band Matrix */
  void zhbmv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *k,
               const std::complex<double> *alpha, const std::complex<double> *a,
               const CPPL_INT *lda, const std::complex<double> *x, const CPPL_INT *incx,
               const std::complex<double> *beta, std::complex<double> *y,
               const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for Hermitian Packed Matrix */
  void zhpmv_( const char *uplo, const CPPL_INT *N,
               const std::complex<double> *alpha, const std::complex<double> *ap,
               const std::complex<double> *x, const CPPL_INT *incx,
               const std::complex<double> *beta, std::complex<double> *y,
               const CPPL_INT *incy );
  
  /* x := A * x for Triangular Matrix */
  void ztrmv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const std::complex<double> *a, const CPPL_INT *lda,
               std::complex<double> *x, const CPPL_INT *incx );
  /* x := A * x for Triangular (Banded Storage) Matrix */
  void ztbmv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const CPPL_INT *k, const std::complex<double> *a,
               const CPPL_INT *lda, std::complex<double> *x, const CPPL_INT *incx );
  /* x := A * x for Triangular (Packed Storage) Matrix */
  void ztpmv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const std::complex<double> *ap,
               std::complex<double> *x, const CPPL_INT *incx );
  
  /* Solve A * x = b for Triangular Matrix */
  void ztrsv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const std::complex<double> *a, const CPPL_INT *lda,
               std::complex<double> *x, const CPPL_INT *incx );
  /* Solve A * x = b for Triangular (Banded Storage) Matrix */
  void ztbsv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const CPPL_INT *k, const std::complex<double> *a,
               const CPPL_INT *lda, std::complex<double> *x, const CPPL_INT *incx );
  /* Solve A * x = b for Triangular (Packed Storage) Matrix */
  void ztpsv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const std::complex<double> *ap,
               std::complex<double> *x, const CPPL_INT *incx );
  
  /* A := alpha * x * y' + A  (A: M by N Matrix) */
  void zgeru_( const CPPL_INT *M, const CPPL_INT *N, const std::complex<double> *alpha,
               const std::complex<double> *x, const CPPL_INT *incx,
               const std::complex<double> *y, const CPPL_INT *incy,
               std::complex<double> *a, const CPPL_INT *lda );
  /* A := alpha * x * conjg(y') + A  (A: M by N Matrix) */
  void zgerc_( const CPPL_INT *M, const CPPL_INT *N, const std::complex<double> *alpha,
               const std::complex<double> *x, const CPPL_INT *incx,
               const std::complex<double> *y, const CPPL_INT *incy,
               std::complex<double> *a, const CPPL_INT *lda );
  /* A := alpha * x * conjg(x') + A  (A: Hermitian N by N Matrix) */
  void zher_( const char *uplo, const CPPL_INT *N, const double *alpha,
              const std::complex<double> *x, const CPPL_INT *incx,
              std::complex<double> *a, const CPPL_INT *lda );
  /* A := alpha * x * conjg(x') + A
     (A: Hermitian N by N Packed Storage Matrix) */
  void zhpr_( const char *uplo, const CPPL_INT *N, const double *alpha,
              const std::complex<double> *x, const CPPL_INT *incx,
              std::complex<double> *ap );
  /* A := alpha * x * conjg(y') + conjg(alpha) * y * conjg(x') + A
     (A: Hermitian N by N Matrix) */
  void zher2_( const char *uplo, const CPPL_INT *N,
               const std::complex<double> *alpha, const std::complex<double> *x,
               const CPPL_INT *incx, const std::complex<double> *y,
               const CPPL_INT *incy, std::complex<double> *a, const CPPL_INT *lda );
  /* A := alpha * x * conjg(y') + conjg(alpha) * y * conjg(x') + A
     (A: Hermitian N by N Packed Matrix) */
  void zhpr2_( const char *uplo, const CPPL_INT *N,
               const std::complex<double> *alpha, const std::complex<double> *x,
               const CPPL_INT *incx, const std::complex<double> *y,
               const CPPL_INT *incy, std::complex<double> *ap );
  
  // Level 3 BLAS
  /* C := alpha * A * B + beta * C  (C: M by N Matrix) */
  void zgemm_( const char *transa, const char *transb, const CPPL_INT *M,
               const CPPL_INT *N, const CPPL_INT *k, const std::complex<double> *alpha,
               const std::complex<double> *a, const CPPL_INT *lda,
               const std::complex<double> *b, const CPPL_INT *ldb,
               const std::complex<double> *beta, std::complex<double> *c,
               const CPPL_INT *ldc );
  /* C := alpha * A * B + beta * C
     (A: Symmetric Matrix,  B, C: M by N Matrices) */
  void zsymm_( const char *side, const char *uplo, const CPPL_INT *M,
               const CPPL_INT *N, const std::complex<double> *alpha,
               const std::complex<double> *a, const CPPL_INT *lda,
               const std::complex<double> *b, const CPPL_INT *ldb,
               const std::complex<double> *beta, std::complex<double> *c,
               const CPPL_INT *ldc );
  /* C := alpha * A * B + beta * C
     (A: Hermitian Matrix,  B, C: M by N Matrices) */
  void zhemm_( const char *side, const char *uplo, const CPPL_INT *M,
               const CPPL_INT *N, const std::complex<double> *alpha,
               const std::complex<double> *a, const CPPL_INT *lda,
               const std::complex<double> *b, const CPPL_INT *ldb,
               const std::complex<double> *beta, std::complex<double> *c,
               const CPPL_INT *ldc );
  /* C:= alpha * A * A' + beta * C
     (A: N by k Matrix,  C: Symmetric Matrix) */
  void zsyrk_( const char *uplo, const char *trans, const CPPL_INT *N,
               const CPPL_INT *k, const std::complex<double> *alpha,
               const std::complex<double> *a, const CPPL_INT *lda,
               const std::complex<double> *beta, std::complex<double> *c,
               const CPPL_INT *ldc );
  /* C:= alpha * A * conjg(A') + beta * C
     (A: N by k Matrix,  C: Hermirian Matrix) */
  void zherk_( const char *uplo, const char *trans, const CPPL_INT *N,
               const CPPL_INT *k, const double *alpha,
               const std::complex<double> *a, const CPPL_INT *lda,
               const double *beta, std::complex<double> *c, const CPPL_INT *ldc );
  /* C := alpha * A * B' + alpha * B * A' + beta * C
     (A, B: N by k Marices,  C: Symmetric Matrix ) */
  void zsyr2k_( const char *uplo, const char *trans, const CPPL_INT *N,
                const CPPL_INT *k, const std::complex<double> *alpha,
                const std::complex<double> *a, const CPPL_INT *lda,
                const std::complex<double> *b, const CPPL_INT *ldb,
                const std::complex<double> *beta, std::complex<double> *c,
                const CPPL_INT *ldc );
  /* C := alpha * A * conjg(B') + conjg(alpha) * B * conjg(A') + beta * C
     (A, B: N by k Marices,  C: Hermitian Matrix ) */
  void zher2k_( const char *uplo, const char *trans, const CPPL_INT *N,
                const CPPL_INT *k, const std::complex<double> *alpha,
                const std::complex<double> *a, const CPPL_INT *lda,
                const std::complex<double> *b, const CPPL_INT *ldb,
                const double *beta, std::complex<double> *c, const CPPL_INT *ldc );
  /* B := alpha * A * B  (A: Triangular Matrix,  B: M by N Matrix) */
  void ztrmm_( const char *side, const char *uplo, const char *transa,
               const char *diag, const CPPL_INT *M, const CPPL_INT *N,
               const std::complex<double> *alpha, const std::complex<double> *a,
               const CPPL_INT *lda, std::complex<double> *b, const CPPL_INT *ldb );
  /* Solve A * X = alpha * B
     (A: Triangular Matrix,  X, B: M by N Matrices) */
  void ztrsm_( const char *side, const char *uplo, const char *transa,
               const char *diag, const CPPL_INT *M, const CPPL_INT *N,
               const std::complex<double> *alpha, const std::complex<double> *a,
               const CPPL_INT *lda, std::complex<double> *b, const CPPL_INT *ldb );
}
