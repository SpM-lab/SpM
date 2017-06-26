extern "C" {
  // Solve Linear Equations  A * x = b
  /* for General Matrix */
  void dgesv_( const CPPL_INT *N, const CPPL_INT *nrhs, double *a, const CPPL_INT *lda,
               CPPL_INT *ipiv, double *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for General Band Matrix */
  void dgbsv_( const CPPL_INT *N, const CPPL_INT *KL, const CPPL_INT *KU,
               const CPPL_INT *nrhs, double *ab, const CPPL_INT *ldab,
               CPPL_INT *ipiv, double *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for Tridiagonal Matrix */
  void dgtsv_( const CPPL_INT *N, const CPPL_INT *nrhs, double *dl, double *d,
               double *du, double *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for Symmetric Positive Definite Matrix */
  void dposv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               double *a, const CPPL_INT *lda, double *b, const CPPL_INT *ldb,
               CPPL_INT *info );
  /* for Symmetric Positive Definite (Packed Storage) Matrix */
  void dppsv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               double *ap, double *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for Symmetric Positive Definite Band Matrix */
  void dpbsv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *kd,
               const CPPL_INT *nrhs, double *ab, const CPPL_INT *ldab,
               double *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for Symmetric Positive Definite Tridiagonal Matrix */
  void dptsv_( const CPPL_INT *N, const CPPL_INT *nrhs, double *d, double *e,
               double *b, const CPPL_INT *ldb, CPPL_INT *info );
  /* for Symmetric Indefinite Matrix */
  void dsysv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               double *a, const CPPL_INT *lda, CPPL_INT *ipiv, double *b,
               const CPPL_INT *ldb, double *work, const CPPL_INT *lwork,
               CPPL_INT *info );
  /* for Symmetric Indefinite (Packed Storage) Matrix */
  void dspsv_( const char *uplo, const CPPL_INT *N, const CPPL_INT *nrhs,
               double *ap, CPPL_INT *ipiv, double *b, const CPPL_INT *ldb,
               CPPL_INT *info );

  // Linear Least Square Problems
  // Solve Overdetermined or Underdetermined Linear Equations
  /* Using Orthogonal Factorization, Assuming Full Rank */
  void dgels_( const char *trans, const CPPL_INT *M, const CPPL_INT *N,
               const CPPL_INT *nrhs, double *a, const CPPL_INT *lda,
               double *b, const CPPL_INT *ldb,
               double *work, const CPPL_INT *lwork, CPPL_INT *info );
  /* Compute Minimum-Norm Solution using Orthogonal Factorization */
  void dgelsy_( const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *nrhs,
                double *a, const CPPL_INT *lda, double *b, const CPPL_INT *ldb,
                CPPL_INT *jpvt, const double *rcond, CPPL_INT *rank,
                double *work, const CPPL_INT *lwork, CPPL_INT *info );
  /* Compute Minimum-Norm Solution using Singular Value Decomposition */
  void dgelss_( const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *nrhs,
                double *a, const CPPL_INT *lda, double *b, const CPPL_INT *ldb,
                double *s, const double *rcond, CPPL_INT *rank,
                double *work, const CPPL_INT *lwork, CPPL_INT *info );
  /* Compute Minimum-Norm Solution using Singular Value Decomposition with Divide and Conquer Algorithm */
  void dgelsd_( const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *nrhs,
                double *a, const CPPL_INT *lda, double *b, const CPPL_INT *ldb,
                double *s, const double *rcond, CPPL_INT *rank,
                double *work, const CPPL_INT *lwork, CPPL_INT *iwork, CPPL_INT *info );
  /* Solve Linear Equality-Constrained Least Squares (LSE) Problem */
  void dgglse_( const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *p, double *a,
                const CPPL_INT *lda, double *b, const CPPL_INT *ldb,
                double *c, double *d, double *x, double *work,
                const CPPL_INT *lwork, CPPL_INT *info );
  /* Solve Gauss-Markov Linear Model (GLM) Problem */
  void dggglm_( const CPPL_INT *N, const CPPL_INT *M, const CPPL_INT *p,
                double *a, const CPPL_INT *lda, double *b, const CPPL_INT *ldb,
                double *d, double *x, double *y,
                double *work, const CPPL_INT *lwork, CPPL_INT *info );

  // Standard Eigenvalue and Singular Value Problems
  // Divide and Conquer routines are more fast, but need more memory.
  /* Eigenvalues/Eigenvectors for General Matrix */
  void dgeev_( const char *jobvl, const char *jobvr, const CPPL_INT *N,
               double *a, const CPPL_INT *lda, double *wr, double *wi,
               double *vl, const CPPL_INT *ldvl, double *vr, const CPPL_INT *ldvr,
               double *work, const CPPL_INT *lwork, CPPL_INT *info );
  /* Eigenvalues/Eigenvectors for Symmetric Matrix */
  void dsyev_( const char *jobz, const char *uplo, const CPPL_INT *N,
               double *a, const CPPL_INT *lda, double *w, double *work,
               const CPPL_INT *lwork, CPPL_INT *info );
  void dsyevd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                double *a, const CPPL_INT *lda, double *w, double *work,
                const CPPL_INT *lwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Eigenvalues/Eigenvectors for Symmetric (Packed Storage) Matrix */
  void dspev_( const char *jobz, const char *uplo, const CPPL_INT *N,
               double *ap, double *w, double *z, const CPPL_INT *ldz,
               double *work, CPPL_INT *info );
  void dspevd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                double *ap, double *w, double *z, const CPPL_INT *ldz,
                double *work, const CPPL_INT *lwork, CPPL_INT *iwork,
                const CPPL_INT *liwork, CPPL_INT *info );  /* Divide and Conqure */
  /* Eigenvalues/Eigenvectors for Symmetric Band Matrix */
  void dsbev_( const char *jobz, const char *uplo, const CPPL_INT *N,
               const CPPL_INT *kd, double *ab, const CPPL_INT *ldab, double *w,
               double *z, const CPPL_INT *ldz, double *work, CPPL_INT *info );
  void dsbevd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                const CPPL_INT *kd, double *ab, const CPPL_INT *ldab, double *w,
                double *z, const CPPL_INT *ldz, double *work,
                const CPPL_INT *lwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Eigenvalues/Eigenvectors for Symmetric Tridiagonal Matrix */
  void dstev_( const char *jobz, const CPPL_INT *N, double *d, double *e,
               double *z, const CPPL_INT *ldz, double *work, CPPL_INT *info );
  void dstevd_( const char *jobz, const CPPL_INT *N, double *d, double *e,
                double *z, const CPPL_INT *ldz, double *work,
                const CPPL_INT *lwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Schur Factorization for General Matrix */
  void dgees_( const char *jobvs, const char *sort,
               bool (*select)( double *, double * ),
               const CPPL_INT *N, double *a, const CPPL_INT *lda, CPPL_INT *sdim,
               double *wr, double *wi, double *vs, const CPPL_INT *ldvs,
               double *work, const CPPL_INT *lwork, bool *bwork,
               CPPL_INT *info );
  /* Singular Value Decomposition for General Matrix */
  void dgesvd_( const char *jobu, const char *jobvt, const CPPL_INT *M,
                const CPPL_INT *N, double *a, const CPPL_INT *lda, double *s,
                double *u, const CPPL_INT *ldu, double *vt, const CPPL_INT *ldvt,
                double *work, const CPPL_INT *lwork, CPPL_INT *info );
  void dgesdd_( const char *jobz, const CPPL_INT *M, const CPPL_INT *N, double *a,
                const CPPL_INT *lda, double *s, double *u, const CPPL_INT *ldu,
                double *vt, const CPPL_INT *ldvt, double *work,
                const CPPL_INT *lwork, CPPL_INT *iwork,
                CPPL_INT *info );  /* Divide and Conqure */

  // Generalized Eigenvalue and Sigular Value Problems
  /* Generalized Eigenvalues/Eigenvectors for General Matrix */
  void dggev_( const char *jobvl, const char *jobvr, const CPPL_INT *N,
               double *a, const CPPL_INT *lda, double *b, const CPPL_INT *ldb,
               double *alphar, double *alphai, double *beta,
               double *vl, const CPPL_INT *ldvl, double *vr, const CPPL_INT *ldvr,
               double *work, const CPPL_INT *lwork, CPPL_INT *info );
  /* Generalized Eigenvalues/Eigenvectors
     for Symmetric-definite Matrix */
  void dsygv_( const CPPL_INT *itype, const char *jobz, const char *uplo,
               const CPPL_INT *N, double *a, const CPPL_INT *lda, double *b,
               const CPPL_INT *ldb, double *w, double *work, const CPPL_INT *lwork,
               CPPL_INT *info );
  void dsygvd_( const CPPL_INT *itype, const char *jobz, const char *uplo,
                const CPPL_INT *N, double *a, const CPPL_INT *lda, double *b,
                const CPPL_INT *ldb, double *w, double *work,
                const CPPL_INT *lwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Generalized Eigenvalues/Eigenvectors
     for Symmetric-definite (Packed Storage) Matrix */
  void dspgv_( const CPPL_INT *itype, const char *jobz, const char *uplo,
               const CPPL_INT *N, double *ap, double *bp, double *w,
               double *z, const CPPL_INT *ldz, double *work, CPPL_INT *info );
  void dspgvd_( const CPPL_INT *itype, const char *jobz, const char *uplo,
                const CPPL_INT *N, double *ap, double *bp, double *w,
                double *z, const CPPL_INT *ldz, double *work,
                const CPPL_INT *lwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Generalized Eigenvalues/Eigenvectors
     Symmetric-definite Band Matrix */
  void dsbgv_( const char *jobz, const char *uplo, const CPPL_INT *N,
               const CPPL_INT *ka, const CPPL_INT *kb, double *ab, const CPPL_INT *ldab,
               double *bb, const CPPL_INT *ldbb, double *w, double *z,
               const CPPL_INT *ldz, double *work, CPPL_INT *info );
  void dsbgvd_( const char *jobz, const char *uplo, const CPPL_INT *N,
                const CPPL_INT *ka, const CPPL_INT *kb, double *ab,
                const CPPL_INT *ldab, double *bb, const CPPL_INT *ldbb, double *w,
                double *z, const CPPL_INT *ldz, double *work,
                const CPPL_INT *lwork, CPPL_INT *iwork, const CPPL_INT *liwork,
                CPPL_INT *info );  /* Divide and Conqure */
  /* Generalized Schur Factrization for General Matrix */
  void dgges_( const char *jobvsl, const char *jobvsr, const char *sort,
               CPPL_INT (*delctg)( double *, double *, double * ),
               const CPPL_INT *N, double *a, const CPPL_INT *lda, double *b,
               const CPPL_INT *ldb, CPPL_INT *sdim, double *alphar, double *alphai,
               double *beta, double *vsl, const CPPL_INT *ldvsl, double *vsr,
               const CPPL_INT *ldvsr, double *work, const CPPL_INT *lwork,
               bool *bwork, CPPL_INT *info );
  /* Generailized Singular Value Decomposition for General Matrix */
  void dggsvd_( const char *jobu, const char *jobv, const char *jobq,
                const CPPL_INT *M, const CPPL_INT *N, const CPPL_INT *p, CPPL_INT *k,
                CPPL_INT *L, double *a, const CPPL_INT *lda, double *b,
                const CPPL_INT *ldb, double *alpha, double *beta,
                double *u, const CPPL_INT *ldu, double *v, const CPPL_INT *ldv,
                double *q, const CPPL_INT *ldq, double *work, CPPL_INT *iwork,
                CPPL_INT *info );
}
