//=============================================================================
/*! solve A*X=Y using dsysv\n
  The argument is dmatrix Y. Y is overwritten and become the solution X.
  A is also overwritten. 
*/
inline CPPL_INT dsymatrix::dsysv(dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=mat.n){
    ERROR_REPORT;
    std::cerr << "These two matrices cannot be solved." << std::endl
              << "Your input was (" << n << "x" << n << ") and (" << mat.n << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  char UPLO('l');
  CPPL_INT NRHS(mat.n), LDA(n), *IPIV(new CPPL_INT[n]), LDB(mat.m), LWORK(-1), INFO(1);
  double *WORK( new double[1] );
  dsysv_(&UPLO, &n, &NRHS, array, &LDA, IPIV, mat.array, &LDB, WORK, &LWORK, &INFO);
  
  INFO=1; 
  LWORK = CPPL_INT(WORK[0]);
  delete [] WORK;
  WORK = new double[LWORK];
  dsysv_(&UPLO, &n, &NRHS, array, &LDA, IPIV, mat.array, &LDB, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] IPIV;
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! solve A*x=y using dsysv\n
  The argument is dcovector y. y is overwritten and become the solution x.
  A is also overwritten.
*/
inline CPPL_INT dsymatrix::dsysv(dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(n!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << n << "x" << n << ") and (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG 
  
  char UPLO('l');
  CPPL_INT NRHS(1), LDA(n), *IPIV(new CPPL_INT[n]), LDB(vec.l), LWORK(-1), INFO(1);
  double *WORK( new double[1] );
  dsysv_(&UPLO, &n, &NRHS, array, &LDA, IPIV, vec.array, &LDB, WORK, &LWORK, &INFO);
  
  INFO=1;
  LWORK = CPPL_INT(WORK[0]);
  delete [] WORK;
  WORK = new double[LWORK];
  dsysv_(&UPLO, &n, &NRHS, array, &LDA, IPIV, vec.array, &LDB, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] IPIV;

  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! calculate eigenvalues and eigenvectors.\n
  All of the arguments need not to be initialized.
  w is overwitten and become eigenvalues.
  This matrix is also overwritten. 
  if jobz=1, this matrix becomes eigenvectors.
*/
inline CPPL_INT dsymatrix::dsyev(std::vector<double>& w, const bool& jobz=0)
{CPPL_VERBOSE_REPORT;
  w.resize(n);
  char JOBZ, UPLO('l');
  if(jobz==0){ JOBZ='n'; } else{ JOBZ='V'; }
  CPPL_INT LDA(n), INFO(1), LWORK(-1);
  double *WORK(new double[1]);
  dsyev_(&JOBZ, &UPLO, &n, array, &LDA, &w[0], WORK, &LWORK, &INFO);
  
  INFO=1;
  LWORK = CPPL_INT(WORK[0]);
  delete [] WORK;
  WORK = new double[LWORK];
  dsyev_(&JOBZ, &UPLO, &n, array, &LDA, &w[0], WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate eigenvalues and eigenvectors.\n
  All of the arguments need not to be initialized.
  w and v are overwitten and become 
  eigenvalues and eigenvectors, respectively.
  This matrix is also overwritten. 
*/
inline CPPL_INT dsymatrix::dsyev(std::vector<double>& w, std::vector<dcovector>& v)
{CPPL_VERBOSE_REPORT;
  w.resize(n);
  v.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    v[i].resize(n);
  }
  
  char JOBZ('V'), UPLO('l');
  CPPL_INT LDA(n), INFO(1), LWORK(-1);
  double *WORK(new double[1]);
  dsyev_(&JOBZ, &UPLO, &n, array, &LDA, &w[0], WORK, &LWORK, &INFO);
  
  INFO=1;
  LWORK = CPPL_INT(WORK[0]);
  delete [] WORK;
  WORK = new double[LWORK];
  dsyev_(&JOBZ, &UPLO, &n, array, &LDA, &w[0], WORK, &LWORK, &INFO);
  delete [] WORK;
  
  //// forming ////
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=0; i<n; i++){
      v[j](i) = array[i+n*j];
    }
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate eigenvalues and eigenvectors.\n
  All of the arguments need not to be initialized.
  w and v are overwitten and become 
  eigenvalues and eigenvectors, respectively.
  This matrix is also overwritten. 
*/
inline CPPL_INT dsymatrix::dsyev(std::vector<double>& w, std::vector<drovector>& v)
{CPPL_VERBOSE_REPORT;
  w.resize(n);
  v.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    v[i].resize(n);
  }
  
  char JOBZ('V'), UPLO('l');
  CPPL_INT LDA(n), INFO(1), LWORK(-1);
  double *WORK(new double[1]);
  dsyev_(&JOBZ, &UPLO, &n, array, &LDA, &w[0], WORK, &LWORK, &INFO);
  
  INFO=1;
  LWORK = CPPL_INT(WORK[0]);
  delete [] WORK;
  WORK = new double[LWORK];
  dsyev_(&JOBZ, &UPLO, &n, array, &LDA, &w[0], WORK, &LWORK, &INFO);
  delete [] WORK;
  
  //// forming ////
  for(CPPL_INT j=0; j<n; j++){
    for(CPPL_INT i=0; i<n; i++){
      v[j](i) = array[i+n*j];
    }
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! calculate generalized eigenvalues\n
  w is overwitten and become generalized eigenvalues.
  This matrix and matB are also overwritten. 
*/
inline CPPL_INT dsymatrix::dsygv(dsymatrix& matB, std::vector<double>& w)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matB.n!=n){
    ERROR_REPORT;
    std::cerr << "The matrix B is not a matrix having the same size as \"this\" matrix." << std::endl
              << "The B matrix is (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  w.resize(n);
  char JOBZ('n'), UPLO('l');
  CPPL_INT ITYPE(1), LDA(n), LDB(n), LWORK(-1), INFO(1);
  double *WORK(new double[1]);
  dsygv_(&ITYPE, &JOBZ, &UPLO, &n, array, &LDA, matB.array, &LDB, &w[0], WORK, &LWORK, &INFO);
  
  INFO=1;
  LWORK = CPPL_INT(WORK[0]);
  delete [] WORK;
  WORK = new double[LWORK];
  dsygv_(&ITYPE, &JOBZ, &UPLO, &n, array, &LDA, matB.array, &LDB, &w[0], WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate generalized eigenvalues\n
  w is overwitten and become generalized eigenvalues.
  This matrix and matB are also overwritten. 
*/
inline CPPL_INT dsymatrix::dsygv(dsymatrix& matB, std::vector<double>& w,
                             std::vector<dcovector>& v)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(matB.n!=n){
    ERROR_REPORT;
    std::cerr << "The matrix B is not a matrix having the same size as \"this\" matrix." << std::endl
              << "The B matrix is (" << matB.n << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  w.resize(n);
  v.resize(n);
  char JOBZ('V'), UPLO('l');
  CPPL_INT ITYPE(1), LDA(n), LDB(n), LWORK(-1), INFO(1);
  double *WORK(new double[1]);
  dsygv_(&ITYPE, &JOBZ, &UPLO, &n, array, &LDA, matB.array, &LDB, &w[0], WORK, &LWORK, &INFO);
  
  INFO=1;
  LWORK = CPPL_INT(WORK[0]);
  delete [] WORK;
  WORK = new double[LWORK];
  dsygv_(&ITYPE, &JOBZ, &UPLO, &n, array, &LDA, matB.array, &LDB, &w[0], WORK, &LWORK, &INFO);
  delete [] WORK;
  
  //// reforming ////
  for(int i=0; i<n; i++){
    v[i].resize(n);
    for(int j=0; j<n; j++){
      v[i](j) =darray[i][j];
    }
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}
