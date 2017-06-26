extern "C" {
  // Solve Linear Equations  A * x = b
  /* for General Matrix */
  void zgesv_( const CPPL_INT *N, const CPPL_INT *nrhs, std::complex<double> *a,
               const CPPL_INT *lda, CPPL_INT *ipiv, std::complex<double> *b,
               const CPPL_INT *ldb, CPPL_INT *info );
  /* for General Band Matrix */
  void zgbsv_( const CPPL_INT *N, const CPPL_INT *KL, const CPPL_INT *KU,
               const CPPL_INT *nrhs, std::complex<double> *ab, const CPPL_INT *ldab,
               CPPL_INT *ipiv, std::complex<double> *b, const CPPL_INT *ldb,
               CPPL_INT *info );
  /* for Tridiagonal Matrix */
  void zgtsv_( const CPPL_INT *N, const CPPL_INT *nrhs, std::complex<double> *dl,
               std::complex<double> *d, std::complex<double> *du,
               std::complex<double> *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for Hermitian Positive Definite Matrix */
  void zposv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               std::complex<double> *a, const CPPL_INT *lda, std::complex<double> *b,
               const CPPL_INT *ldb, CPPL_INT *info );
  /* for Hermitian Positive Definite (Packed Storage) Matrix */
  void zppsv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               std::complex<double> *ap, std::complex<double> *b, const CPPL_INT *ldb,
               CPPL_INT *info );
  /* for Hermitian Positive Definite Band Matrix */
  void zpbsv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *kd,
               const CPPL_INT *nrhs, std::complex<double> *ab, const CPPL_INT *ldab,
               std::complex<double> *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for Hermitian Positive Definite Tridiagonal Matrix */
  void zptsv_( const CPPL_INT *N, const CPPL_INT *nrhs, double *d,
               std::complex<double> *e, std::complex<double> *b, const CPPL_INT *ldb,
               CPPL_INT *info );
  /* for Hermitian Indefinite Matrix */
  void zhesv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               std::complex<double> *a, const CPPL_INT *lda, CPPL_INT *ipiv,
               std::complex<double> *b, const CPPL_INT *ldb, std::complex<double> *work,
               const CPPL_INT *lwork, CPPL_INT *info );
  /* for Symmetric Indefinite Matrix */
  void zsysv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               std::complex<double> *a, const CPPL_INT *lda, CPPL_INT *ipiv,
               std::complex<double> *b, const CPPL_INT *ldb, std::complex<double> *work,
               const CPPL_INT *lwork, CPPL_INT *info );
  /* for Hermitian Indefinite (Packed Storage) Matrix */
  void zhpsv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               std::complex<double> *ap, CPPL_INT *ipiv, std::complex<double> *b,
               const CPPL_INT *ldb, CPPL_INT *info );
  /* for Symmetric Indefinite (Packed Storage) Matrix */
  void zspsv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               std::complex<double> *ap, CPPL_INT *ipiv, std::complex<double> *b,
               const CPPL_INT *ldb, CPPL_INT *info );

  // Linear Least Square Problems
  // Solve Overdetermined or Underdetermined Linear Equations
  /* Using Orthogonal Factorization, Assuming Full Rank */
  void zgels_( const char *trans, const CPPL_INT *M, const CPPL_INT *N,
               const CPPL_INT *nrhs, std::complex<double> *a, const CPPL_INT *lda,
               std::complex<double> *b, const CPPL_INT *ldb,
               std::complex<double> *work, const CPPL_INT *lwork, CPPL_INT *info );
  /* Compute Minimum-Norm Solution using Orthogonal Factorization */
  void zgelsy_( const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *nrhs,
                std::complex<double> *a, const CPPL_INT *lda, std::complex<double> *b,
                const CPPL_INT *ldb, CPPL_INT *jpvt, const double *rcond,
                CPPL_INT *rank, std::complex<double> *work, const CPPL_INT *lwork,
                double *rwork, CPPL_INT *info );
  /* Compulte Minimum-Norm Solution using Singular Value Decomposition */
  void zgelss_( const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *nrhs,
                std::complex<double> *a, const CPPL_INT *lda, std::complex<double> *b,
                const CPPL_INT *ldb, double *s, const double *rcond,
                CPPL_INT *rank, std::complex<double> *work, const CPPL_INT *lwork,
                double *rwork, CPPL_INT *info );
  /* Solve Linear Equality-Constrained Least Squares (LSE) Problem */
  void zgglse_( const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *p,
                std::complex<double> *a, const CPPL_INT *lda, std::complex<double> *b,
                const CPPL_INT *ldb, std::complex<double> *c, std::complex<double> *d,
                std::complex<double> *x, std::complex<double> *work,
                const CPPL_INT *lwork, CPPL_INT *info );
  /* Solve Gauss-Markov Linear Model (GLM) Problem */
  void zggglm_( const CPPL_INT *N, const CPPL_INT *M, const CPPL_INT *p,
                std::complex<double> *a, const CPPL_INT *lda, std::complex<double> *b,
                const CPPL_INT *ldb, std::complex<double> *d, std::complex<double> *x,
                std::complex<double> *y, std::complex<double> *work,
                const CPPL_INT *lwork, CPPL_INT *info );

  // Standard Eigenvalue and Singular Value Problems
  // Divide and Conquer routines are more fast, but need more memory.
  /* Eigenvalues/Eigenvectors for General Matrix */
  void zgeev_( const char *jobvl, const char *jobvr, const CPPL_INT *N,
               std::complex<double> *a, const CPPL_INT *lda, std::complex<double> *w,
               std::complex<double> *vl, const CPPL_INT *ldvl, std::complex<double> *vr,
               const CPPL_INT *ldvr, std::complex<double> *work, const CPPL_INT *lwork,
               double *rwork, CPPL_INT *info );
  /* Eigenvalues/Eigenvectors for Hermitian Matrix */
  void zheev_( const char *jobz, const char *uplo, const CPPL_INT *N,
               std::complex<double> *a, const CPPL_INT *lda, double *w,
               std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
               CPPL_INT *info );
  void zheevd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                std::complex<double> *a, const CPPL_INT *lda, double *w,
                std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
                const CPPL_INT *lrwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Eigenvalues/Eigenvectors for Hermitian (Packed Storage) Matrix */
  void zhpev_( const char *jobz, const char *uplo, const CPPL_INT *N,
               std::complex<double> *ap, double *w, std::complex<double> *z,
               const CPPL_INT *ldz, std::complex<double> *work, double *rwork,
               CPPL_INT *info );
  void zhpevd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                std::complex<double> *ap, double *w, std::complex<double> *z,
                const CPPL_INT *ldz, std::complex<double> *work, const CPPL_INT *lwork,
                double *rwork, const CPPL_INT *lrwork, CPPL_INT *iwork,
                const CPPL_INT *liwork, CPPL_INT *info );  /* Divide and Conqure */
  /* Eigenvalues/Eigenvectors for Hermitian Band Matrix */
  void zhbev_( const char *jobz, const char *uplo, const CPPL_INT *N,
               const CPPL_INT *kd, std::complex<double> *ab, const CPPL_INT *ldab,
               double *w, std::complex<double> *z, const CPPL_INT *ldz,
               std::complex<double> *work, double *rwork, CPPL_INT *info );
  void zhbevd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                const CPPL_INT *kd, std::complex<double> *ab, const CPPL_INT *ldab,
                double *w, std::complex<double> *z, const CPPL_INT *ldz,
                std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
                const CPPL_INT *lrwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Schur Factorization for General Matrix */
  void zgees_( const char *jobvs, const char *sort,
               bool (*select)( std::complex<double> *, std::complex<double> * ),
               const CPPL_INT *N, std::complex<double> *a, const CPPL_INT *lda,
               CPPL_INT *sdim, std::complex<double> *w, std::complex<double> *vs,
               const CPPL_INT *ldvs, std::complex<double> *work, const CPPL_INT *lwork,
               double *rwork, bool *bwork, CPPL_INT *info );
  /* Singular Value Decomposition for General Matrix */
  void zgesvd_( const char *jobu, const char *jobvt, const CPPL_INT *M,
                const CPPL_INT *N, std::complex<double> *a, const CPPL_INT *lda,
                double *s, std::complex<double> *u, const CPPL_INT *ldu,
                std::complex<double> *vt, const CPPL_INT *ldvt,
                std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
                CPPL_INT *info );
  void zgesdd_( const char *jobz, const CPPL_INT *M, const CPPL_INT *N,
                std::complex<double> *a, const CPPL_INT *lda, double *s,
                std::complex<double> *u, const CPPL_INT *ldu, std::complex<double> *vt,
                const CPPL_INT *ldvt, std::complex<double> *work, const CPPL_INT *lwork,
                double *rwork, CPPL_INT *iwork,
                CPPL_INT *info );  /* Divide and Conqure */

  // Generalized Eigenvalue and Sigular Value Problems
  /* Generalized Eigenvalues/Eigenvectors for General Matrix */
  void zggev_( const char *jobvl, const char *jobvr, const CPPL_INT *N,
               std::complex<double> *a, const CPPL_INT *lda, std::complex<double> *b,
               const CPPL_INT *ldb, std::complex<double> *alpha,
               std::complex<double> *beta, std::complex<double> *vl,
               const CPPL_INT *ldvl, std::complex<double> *vr, const CPPL_INT *ldvr,
               std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
               CPPL_INT *info );
  /* Generalized Eigenvalues/Eigenvectors
     for Hermitian-definite Matrix */
  void zhegv_( const CPPL_INT *itype, const char *jobz, const char *uplo,
               const CPPL_INT *N, std::complex<double> *a, const CPPL_INT *lda,
               std::complex<double> *b, const CPPL_INT *ldb, double *w,
               std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
               CPPL_INT *info );
  void zhegvd_( const CPPL_INT *itype, const char *jobz, const char *uplo,
                const CPPL_INT *N, std::complex<double> *a, const CPPL_INT *lda,
                std::complex<double> *b, const CPPL_INT *ldb, double *w,
                std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
                const CPPL_INT *lrwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Generalized Eigenvalues/Eigenvectors
     for Hermitian-definite (Packed Storage) Matrix */
  void zhpgv_( const CPPL_INT *itype, const char *jobz, const char *uplo,
               const CPPL_INT *N, std::complex<double> *ap, std::complex<double> *bp,
               double *w, std::complex<double> *z, const CPPL_INT *ldz,
               std::complex<double> *work, double *rwork, CPPL_INT *info );
  void zhpgvd_( const CPPL_INT *itype, const char *jobz, const char *uplo,
                const CPPL_INT *N, std::complex<double> *ap, std::complex<double> *bp,
                double *w, std::complex<double> *z, const CPPL_INT *ldz,
                std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
                const CPPL_INT *lrwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Generalized Eigenvalues/Eigenvectors
     for Hermitian-definite Band Matrix */
  void zhbgv_( const char *jobz, const char *uplo, const CPPL_INT *N,
               const CPPL_INT *ka, const CPPL_INT *kb, std::complex<double> *ab,
               const CPPL_INT *ldab, std::complex<double> *bb, const CPPL_INT *ldbb,
               double *w, std::complex<double> *z, const CPPL_INT *ldz,
               std::complex<double> *work, double *rwork, CPPL_INT *info );
  void zhbgvd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                const CPPL_INT *ka, const CPPL_INT *kb, std::complex<double> *ab,
                const CPPL_INT *ldab, std::complex<double> *bb, const CPPL_INT *ldbb,
                double *w, std::complex<double> *z, const CPPL_INT *ldz,
                std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
                const CPPL_INT *lrwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Generalized Schur Factrization for General Matrix */
  void zgges_( const char *jobvsl, const char *jobvsr, const char *sort,
               bool (*delctg)( std::complex<double> *, std::complex<double> * ),
               const CPPL_INT *N, std::complex<double> *a, const CPPL_INT *lda,
               std::complex<double> *b, const CPPL_INT *ldb, CPPL_INT *sdim,
               std::complex<double> *alpha, std::complex<double> *beta,
               std::complex<double> *vsl, const CPPL_INT *ldvsl,
               std::complex<double> *vsr, const CPPL_INT *ldvsr,
               std::complex<double> *work, const CPPL_INT *lwork, double *rwork,
               bool *bwork, CPPL_INT *info );
  /* Generailized Singular Value Decomposition for General Matrix */
  void zggsvd_( const char *jobu, const char *jobv, const char *jobq,
                const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *p, CPPL_INT *k,
                CPPL_INT *L, std::complex<double> *a, const CPPL_INT *lda,
                std::complex<double> *b, const CPPL_INT *ldb, double *alpha,
                double *beta, std::complex<double> *u, const CPPL_INT *ldu,
                std::complex<double> *v, const CPPL_INT *ldv, std::complex<double> *q,
                const CPPL_INT *ldq, std::complex<double> *work, double *rwork,
                CPPL_INT *iwork, CPPL_INT *info );
}
