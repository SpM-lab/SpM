extern "C" {
  // Level 1 BLAS
  /* Vector rotation: x := c*x[i] + s*y[i], y[i] := - s*x[i] + c*y[i] */
  void drot_( const CPPL_INT *N, double *x, const CPPL_INT *incx, double *y,
              const CPPL_INT *incy, const double *c, const double *s );
  /* x <--> y */
  void dswap_( const CPPL_INT *N, double *x, const CPPL_INT *incx, double *y,
               const CPPL_INT *incy );
  /* x := alpha * x */
  void dscal_( const CPPL_INT *N, const double *alpha, double *x,
               const CPPL_INT *incx );
  /* y := x */
  void dcopy_( const CPPL_INT *N, const double *x, const CPPL_INT *incx,
               double *y, const CPPL_INT *incy );
  /* y := alpha * x + y */
  void daxpy_( const CPPL_INT *N, const double *alpha, const double *x,
               const CPPL_INT *incx, double *y, const CPPL_INT *incy );
  /* return x^T y */
  double ddot_( const CPPL_INT *N, const double *x, const CPPL_INT *incx,
                const double *y, const CPPL_INT *incy );
  /* return N2 norm */
  double dnrm2_( const CPPL_INT *N, const double *x, const CPPL_INT *incx );
  /* return sum of abs(x[i]) */
  double dasum_( const CPPL_INT *N, const double *x, const CPPL_INT *incx );
  /* return the index of element having the largest absolute value */
  CPPL_INT idamax_( const CPPL_INT *N, const double *x, const CPPL_INT *incx );

  // Level 2 BLAS
  /* y := alpha * A * x + beta * y for General M by N Matrix */
  void dgemv_( const char *trans, const CPPL_INT *M, const CPPL_INT *N,
               const double *alpha, const double *a, const CPPL_INT *lda,
               const double *x, const CPPL_INT *incx, const double *beta,
               double *y, const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for General M by N Band Matrix */
  void dgbmv_( const char *trans, const CPPL_INT *M, const CPPL_INT *N,
               const CPPL_INT *KL, const CPPL_INT *KU, const double *alpha,
               const double *a, const CPPL_INT *lda, const double *x,
               const CPPL_INT *incx, const double *beta, double *y,
               const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for Symmetric Matrix */
  void dsymv_( const char *uplo, const CPPL_INT *N, const double *alpha,
               const double *a, const CPPL_INT *lda, const double *x,
               const CPPL_INT *incx, const double *beta, double *y,
               const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for Symmetric Band Matrix */
  void dsbmv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *k,
               const double *alpha, const double *a, const CPPL_INT *lda,
               const double *x, const CPPL_INT *incx, const double *beta,
               double *y, const CPPL_INT *incy );
  /* y := alpha * A * x + beta * y for Symmetric Packed Matrix */
  void dspmv_( const char *uplo, const CPPL_INT *N, const double *alpha,
               const double *ap, const double *x, const CPPL_INT *incx,
               const double *beta, double *y, const CPPL_INT *incy );

  /* x := A * x for Triangular Matrix */
  void dtrmv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const double *a, const CPPL_INT *lda,
               double *x, const CPPL_INT *incx );
  /* x := A * x for Triangular (Banded Storage) Matrix */
  void dtbmv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const CPPL_INT *k, const double *a,
               const CPPL_INT *lda, double *x, const CPPL_INT *incx );
  /* x := A * x for Triangular (Packed Storage) Matrix */
  void dtpmv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const double *ap, double *x,
               const CPPL_INT *incx );

  /* Solve A * x = b for Triangular Matrix */
  void dtrsv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const double *a, const CPPL_INT *lda, double *x,
               const CPPL_INT *incx );
  /* Solve A * x = b for Triangular (Banded Storage) Matrix */
  void dtbsv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const CPPL_INT *k, const double *a,
               const CPPL_INT *lda, double *x, const CPPL_INT *incx );
  /* Solve A * x = b for Triangular (Packed Storage) Matrix */
  void dtpsv_( const char *uplo, const char *trans, const char *diag,
               const CPPL_INT *N, const double *ap, double *x,
               const CPPL_INT *incx );

  /* A := alpha * x * y' + A  (A: M by N Matrix) */
  void dger_( const CPPL_INT *M, const CPPL_INT *N, const double *alpha,
              const double *x, const CPPL_INT *incx, const double *y,
              const CPPL_INT *incy, double *a, const CPPL_INT *lda );
  /* A := alpha * x * x' + A  (A: Symmetric N by N Matrix) */
  void dsyr_( const char *uplo, const CPPL_INT *N, const double *alpha,
              const double *x, const CPPL_INT *incx, double *a,
              const CPPL_INT *lda );
  /* A := alpha * x * x' + A
     (A: Symmetric N by N Packed Storage Matrix) */
  void dspr_( const char *uplo, const CPPL_INT *N, const double *alpha,
              const double *x, const CPPL_INT *incx, double *ap );
  /* A := alpha * x * y' + alpha * y * x' + A
     (A: Symmetric N by N Matrix) */
  void dsyr2_( const char *uplo, const CPPL_INT *N, const double *alpha,
               const double *x, const CPPL_INT *incx, const double *y,
               const CPPL_INT *incy, double *a, const CPPL_INT *lda );
  /* A := alpha * x * y' + alpha * y * x' + A
     (A: Symmetric N by N Packed Storage Matrix) */
  void dspr2_( const char *uplo, const CPPL_INT *N, const double *alpha,
               const double *x, const CPPL_INT *incx, const double *y,
               const CPPL_INT *incy, double *ap );

  // Level 3 BLAS
  /* C := alpha * A * B + beta * C  (C: M by N Matrix) */
  void dgemm_( const char *transa, const char *transb, const CPPL_INT *M,
               const CPPL_INT *N, const CPPL_INT *k, const double *alpha,
               const double *a, const CPPL_INT *lda, const double *b,
               const CPPL_INT *ldb, const double *beta, double *c,
               const CPPL_INT *ldc );
  /* C := alpha * A * B + beta * C
     (A: Symmetric Matrix,  B, C: M by N Matrices) */
  void dsymm_( const char *side, const char *uplo, const CPPL_INT *M,
               const CPPL_INT *N, const double *alpha, const double *a,
               const CPPL_INT *lda, const double *b, const CPPL_INT *ldb,
               const double *beta, double *c, const CPPL_INT *ldc );
  /* C:= alpha * A * A' + beta * C
     (A: N by k Matrix,  C: Symmetric Matrix) */
  void dsyrk_( const char *uplo, const char *trans, const CPPL_INT *N,
               const CPPL_INT *k, const double *alpha, const double *a,
               const CPPL_INT *lda, const double *beta, double *c,
               const CPPL_INT *ldc );
  /* C := alpha * A * B' + alpha * B * A' + beta * C
     (A, B: N by k Marices,  C: Symmetric Matrix ) */
  void dsyr2k_( const char *uplo, const char *trans, const CPPL_INT *N,
                const CPPL_INT *k, const double *alpha, const double *a,
                const CPPL_INT *lda, const double *b, const CPPL_INT *ldb,
                const double *beta, double *c, const CPPL_INT *ldc );
  /* B := alpha * A * B  (A: Triangular Matrix,  B: M by N Matrix) */
  void dtrmm_( const char *side, const char *uplo, const char *transa,
               const char *diag, const CPPL_INT *M, const CPPL_INT *N,
               const double *alpha, const double *a, const CPPL_INT *lda,
               double *b, const CPPL_INT *ldb );
  /* Solve A * X = alpha * B
     (A: Triangular Matrix,  X, B: M by N Matrices) */
  void dtrsm_( const char *side, const char *uplo, const char *transa,
               const char *diag, const CPPL_INT *M, const CPPL_INT *N,
               const double *alpha, const double *a, const CPPL_INT *lda,
               double *b, const CPPL_INT *ldb );
}
