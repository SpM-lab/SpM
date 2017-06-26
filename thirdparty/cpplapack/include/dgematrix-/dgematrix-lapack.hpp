//=============================================================================
/*! solve A*X=Y using dgesv\n
  The argument is dgematrix Y. Y is overwritten and become the solution X.
  A is also overwritten and become P*l*U.
*/
inline CPPL_INT dgematrix::dgesv(dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n || m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrices cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG 
  
  CPPL_INT NRHS(mat.n), LDA(n), *IPIV(new CPPL_INT[n]), LDB(mat.m), INFO(1);
  dgesv_(&n, &NRHS, array, &LDA, IPIV, mat.array, &LDB, &INFO);
  delete [] IPIV;
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! solve A*x=y using dgesv\n
  The argument is dcovector y. y is overwritten and become the solution x.
  A is also overwritten and become P*l*U.
*/
inline CPPL_INT dgematrix::dgesv(dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n || m!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG 
  
  CPPL_INT NRHS(1), LDA(n), *IPIV(new CPPL_INT[n]), LDB(vec.l), INFO(1);
  dgesv_(&n, &NRHS, array, &LDA, IPIV, vec.array, &LDB, &INFO);
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
/*! solve overdetermined or underdetermined A*X=Y using dgels\n
 */
inline CPPL_INT dgematrix::dgels(dgematrix& mat)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrices cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG    
  
  if(m<n){ //underdetermined
    dgematrix tmp(n,mat.n);
    for(CPPL_INT i=0; i<mat.m; i++){
      for(CPPL_INT j=0; j<mat.n; j++){
        tmp(i,j) =mat(i,j);
      }
    }
    mat.clear();
    swap(mat,tmp);
  }
  
  char TRANS('n');
  CPPL_INT NRHS(mat.n), LDA(m), LDB(mat.m), LWORK(std::min(m,n)+std::max(std::min(m,n),NRHS)), INFO(1);
  double *WORK(new double[LWORK]);
  dgels_(&TRANS, &m, &n, &NRHS, array, &LDA, mat.array, &LDB, WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(m>n){ //overdetermined
    dgematrix tmp(n,mat.n);
    for(CPPL_INT i=0; i<tmp.m; i++){
      for(CPPL_INT j=0; j<tmp.n; j++){
        tmp(i,j) =mat(i,j);
      }
    }
    mat.clear();
    swap(mat,tmp);
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! solve overdetermined or underdetermined A*x=y using dgels\n
 */
inline CPPL_INT dgematrix::dgels(dcovector& vec)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG    
  
  if(m<n){ //underdetermined
    dcovector tmp(n);
    for(CPPL_INT i=0; i<vec.l; i++){
      tmp(i)=vec(i);
    }
    vec.clear();
    swap(vec,tmp);
  }
  
  char TRANS('n');
  CPPL_INT NRHS(1), LDA(m), LDB(vec.l), LWORK(std::min(m,n)+std::max(std::min(m,n),NRHS)), INFO(1);
  double *WORK(new double[LWORK]);
  dgels_(&TRANS, &m, &n, &NRHS, array, &LDA, vec.array, &LDB, WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(m>n){ //overdetermined
    dcovector tmp(n);
    for(CPPL_INT i=0; i<tmp.l; i++){
      tmp(i)=vec(i);
    }
    vec.clear();
    swap(vec,tmp);
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! solve overdetermined or underdetermined A*X=Y using dgels
  with the sum of residual squares output\n
  The residual is set as the columnwise sum of residual squares 
  for overdetermined problems
  while it is always zero for underdetermined problems.
*/
inline CPPL_INT dgematrix::dgels(dgematrix& mat, drovector& residual)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=mat.m){
    ERROR_REPORT;
    std::cerr << "These two matrices cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << mat.m << "x" << mat.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  residual.resize(mat.n);
  residual.zero();
  
  if(m<n){ //underdetermined
    dgematrix tmp(n,mat.n);
    for(CPPL_INT i=0; i<mat.m; i++){
      for(CPPL_INT j=0; j<mat.n; j++){
        tmp(i,j) =mat(i,j);
      }
    }
    mat.clear();
    swap(mat,tmp);
  }
  
  char TRANS('n');
  CPPL_INT NRHS(mat.n), LDA(m), LDB(mat.m), LWORK(std::min(m,n)+std::max(std::min(m,n),NRHS)), INFO(1);
  double *WORK(new double[LWORK]);
  dgels_(&TRANS, &m, &n, &NRHS, array, &LDA, mat.array, &LDB, WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(m>n){ //overdetermined
    for(CPPL_INT i=0; i<residual.l; i++){
      for(CPPL_INT j=0; j<m-n; j++){
        residual(i) += std::pow(mat(n+j,i), 2.0);
      }
    }
    
    dgematrix tmp(n,mat.n);
    for(CPPL_INT i=0; i<tmp.m; i++){
      for(CPPL_INT j=0; j<tmp.n; j++){
        tmp(i,j) =mat(i,j);
      }
    }
    mat.clear();
    swap(mat,tmp);
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! solve overdetermined or underdetermined A*x=y using dgels
  with the sum of residual squares output\n
  The residual is set as the sum of residual squares 
  for overdetermined problems
  while it is always zero for underdetermined problems.
*/
inline CPPL_INT dgematrix::dgels(dcovector& vec, double& residual)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=vec.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << vec.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG    
  
  residual=0.0;
  
  if(m<n){ //underdetermined
    dcovector tmp(n);
    for(CPPL_INT i=0; i<vec.l; i++){
      tmp(i)=vec(i);
    }
    vec.clear();
    swap(vec,tmp);
  }
  
  char TRANS('n');
  CPPL_INT NRHS(1), LDA(m), LDB(vec.l), LWORK(std::min(m,n)+std::max(std::min(m,n),NRHS)), INFO(1);
  double *WORK(new double[LWORK]);
  dgels_(&TRANS, &m, &n, &NRHS, array, &LDA, vec.array, &LDB, WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(m>n){ //overdetermined
    for(CPPL_INT i=0; i<m-n; i++){
      residual+=std::pow(vec(n+i),2.0);
    }
    
    dcovector tmp(n);
    for(CPPL_INT i=0; i<tmp.l; i++){
      tmp(i)=vec(i);
    }
    vec.clear();
    swap(vec,tmp);
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
/*! calculate the least-squares-least-norm solution for overdetermined or 
  underdetermined A*x=y using dgelss\n
*/
inline CPPL_INT dgematrix::dgelss
(
 dcovector& B,
 dcovector& S,
 CPPL_INT& RANK,
 const double RCOND =-1.
 )
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=B.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << B.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG    
  
  if(m<n){ //underdetermined
    dcovector tmp(n);
    for(CPPL_INT i=0; i<B.l; i++){
      tmp(i)=B(i);
    }
    B.clear();
    swap(B,tmp);
  }
  
  S.resize(std::min(m,n));
  
  CPPL_INT NRHS(1), LDA(m), LDB(B.l), LWORK(3*std::min(m,n)+std::max(std::max(2*std::min(m,n),std::max(m,n)), NRHS)), INFO(1);
  double *WORK(new double[LWORK]);
  dgelss_(&m, &n, &NRHS, array, &LDA, B.array, &LDB, S.array, &RCOND, &RANK, WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate the least-squares-least-norm solution for overdetermined or 
  underdetermined A*x=y using dgelss\n
*/
inline CPPL_INT dgematrix::dgelss
(
 dgematrix& B,
 dcovector& S,
 CPPL_INT& RANK,
 const double RCOND =-1.
 )
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=B.m){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << B.m << "x" << B.n  << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG    
  
  if(m<n){ //underdetermined
    dgematrix tmp(n,B.n);
    for(CPPL_INT i=0; i<B.m; i++){
      for(CPPL_INT j=0; j<B.n; j++){
        tmp(i,j)=B(i,j);
      }
    }
    B.clear();
    swap(B,tmp);
  }
  
  S.resize(std::min(m,n));
  
  CPPL_INT NRHS(B.n), LDA(m), LDB(B.m), LWORK(3*std::min(m,n)+std::max(std::max(2*std::min(m,n),std::max(m,n)), NRHS)), INFO(1);
  double *WORK(new double[LWORK]);
  dgelss_(&m, &n, &NRHS, array, &LDA, B.array, &LDB, S.array, &RCOND, &RANK, WORK, &LWORK, &INFO);
  delete [] WORK;

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
/*! calculate the least-squares-least-norm solution for overdetermined or underdetermined A*x=y using dgelsd\n
*/
inline CPPL_INT dgematrix::dgelsd
(
 dcovector& B,
 dcovector& S,
 CPPL_INT& RANK,
 const double RCOND =-1.
 )
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=B.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << B.l << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG    
  
  //////// resize ////////
  if(m<n){ //underdetermined
    dcovector tmp(n);
    for(CPPL_INT i=0; i<B.l; i++){
      tmp(i)=B(i);
    }
    B.clear();
    swap(B,tmp);
  }
  
  S.resize(std::min(m,n));
  
  //////// prerun ////////
  CPPL_INT M =m;
  CPPL_INT N =n;
  CPPL_INT NRHS =1;
  CPPL_INT LDA =m;
  CPPL_INT LDB =B.l;
  std::vector<double> WORK(1);
  CPPL_INT LWORK =-1;
  CPPL_INT MINMN =std::min(M,N);
  CPPL_INT SMLSIZ =25;
  CPPL_INT NLVL =CPPL_INT( std::log(double(MINMN)/double(SMLSIZ+1))/std::log(2.) );
  CPPL_INT LIWORK =3*MINMN*NLVL+11*MINMN;
  std::vector<CPPL_INT> IWORK(LIWORK);
  CPPL_INT INFO =1;
  dgelsd_(&M, &N, &NRHS, array, &LDA, B.array, &LDB, S.array, &RCOND, &RANK, &WORK[0], &LWORK, &IWORK[0], &INFO);
  LWORK =CPPL_INT(WORK[0]);
  WORK.resize(LWORK);
  
  //////// run ////////
  INFO =1;
  RANK =1;
  dgelsd_(&M, &N, &NRHS, array, &LDA, B.array, &LDB, S.array, &RCOND, &RANK, &WORK[0], &LWORK, &IWORK[0], &INFO);
  
  //////// info ////////
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
/*! calculate eigenvalues\n
  All of the arguments need not to be initialized. 
  wr and wi are overwitten and become 
  real and imaginary part of eigenvalues, respectively. 
  This matrix is also overwritten. 
*/
inline CPPL_INT dgematrix::dgeev(std::vector<double>& wr, std::vector<double>& wi)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n){
    ERROR_REPORT;
    std::cerr << "This matrix is not a square matrix." << std::endl
              << "This matrix is (" << m << "x" << n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  wr.resize(n);
  wi.resize(n);
  char JOBVL('n'), JOBVR('n');
  CPPL_INT LDA(n), LDVL(1), LDVR(1), LWORK(3*n), INFO(1);
  double *VL(NULL), *VR(NULL), *WORK(new double[LWORK]);
  dgeev_(&JOBVL, &JOBVR, &n, array, &LDA, &wr[0], &wi[0], VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] VL;
  delete [] VR;
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate eigenvalues\n
  All of the arguments need not to be initialized. 
  w are overwitten and become eigenvalues, respectively.
  This matrix is also overwritten. 
*/
inline CPPL_INT dgematrix::dgeev(zcovector& w)
{CPPL_VERBOSE_REPORT;
  //// call dgeev ////
  std::vector<double> wr, wi;
  CPPL_INT INFO =dgeev(wr,wi);
  
  //// assign ////
  w.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    w(i) =comple(wr[i],wi[i]);
  }
  
  return INFO;
}

//=============================================================================
/*! calculate right eigenvalues and right eigenvectors\n
  All of the arguments need not to be initialized. 
  wr, wi, vrr, vri are overwitten and become 
  real and imaginary part of right eigenvalues and right eigenvectors, 
  respectively. 
  This matrix is also overwritten. 
*/
inline CPPL_INT dgematrix::dgeev
(
 std::vector<double>& wr,
 std::vector<double>& wi, 
 std::vector<dcovector>& vrr,
 std::vector<dcovector>& vri
 )
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n){
    ERROR_REPORT;
    std::cerr << "This matrix is not a square matrix." << std::endl
              << "This matrix is (" << m << "x" << n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  wr.resize(n);
  wi.resize(n);
  vrr.resize(n);
  vri.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    vrr[i].resize(n);
    vri[i].resize(n);
  }
  
  dgematrix VR(n,n);
  char JOBVL('n'), JOBVR('V');
  CPPL_INT LDA(n), LDVL(1), LDVR(n), LWORK(4*n), INFO(1);
  double *VL(NULL), *WORK(new double[LWORK]);
  dgeev_(&JOBVL, &JOBVR, &n, array, &LDA, &wr[0], &wi[0], VL, &LDVL, VR.array, &LDVR, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] VL;
  
  //// forming ////
  for(CPPL_INT j=0; j<n; j++){
    if(fabs(wi[j])<DBL_MIN){
      for(CPPL_INT i=0; i<n; i++){
        vrr[j](i) = VR(i,j);
        vri[j](i) = 0.0;
      }
    }
    else{
      for(CPPL_INT i=0; i<n; i++){
        vrr[j](i)   = VR(i,j);
        vri[j](i)   = VR(i,j+1);
        vrr[j+1](i) = VR(i,j);
        vri[j+1](i) =-VR(i,j+1);
      }
      j++;
    }
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate left eigenvalues and left eigenvectors\n
  All of the arguments need not to be initialized. 
  wr, wi, vrr, vri are overwitten and become 
  real and imaginary part of left eigenvalues and left eigenvectors, 
  respectively. 
  This matrix is also overwritten. 
*/
inline CPPL_INT dgematrix::dgeev(std::vector<double>& wr, std::vector<double>& wi, 
                             std::vector<drovector>& vlr, 
                             std::vector<drovector>& vli)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n){
    ERROR_REPORT;
    std::cerr << "This matrix is not a square matrix." << std::endl
              << "This matrix is (" << m << "x" << n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  wr.resize(n);
  wi.resize(n);
  vlr.resize(n);
  vli.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    vlr[i].resize(n);
    vli[i].resize(n);
  }
  
  dgematrix VL(n,n);
  char JOBVL('V'), JOBVR('n');
  CPPL_INT LDA(n), LDVL(n), LDVR(1), LWORK(4*n), INFO(1);
  double *VR(NULL), *WORK(new double[LWORK]);
  dgeev_(&JOBVL, &JOBVR, &n, array, &LDA, &wr[0], &wi[0], VL.array, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] VR;

  //// forming ////
  for(CPPL_INT j=0; j<n; j++){
    if(fabs(wi[j])<DBL_MIN){
      for(CPPL_INT i=0; i<n; i++){
        vlr[j](i) = VL(i,j);
        vli[j](i) = 0.0;
      }
    }
    else{
      for(CPPL_INT i=0; i<n; i++){
        vlr[j](i)   = VL(i,j);
        vli[j](i)   =-VL(i,j+1);
        vlr[j+1](i) = VL(i,j);
        vli[j+1](i) = VL(i,j+1);
      }
      j++;
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
  All of the arguments don't need to be initialized. 
  wr and wi are overwitten and become 
  real and imaginary part of generalized eigenvalues, respectively. 
  This matrix and matB are also overwritten. 
*/
inline CPPL_INT dgematrix::dggev(dgematrix& matB,
                             std::vector<double>& wr, std::vector<double>& wi)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n){
    ERROR_REPORT;
    std::cerr << "This matrix is not a square matrix." << std::endl
              << "This matrix is (" << m << "x" << n << ")." << std::endl;
    exit(1);
  }
  if(matB.m!=n || matB.n!=n){
    ERROR_REPORT;
    std::cerr << "The matrix B is not a square matrix having the same size as \"this\" matrix." << std::endl
              << "The B matrix is (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  wr.resize(n);
  wi.resize(n);
  char JOBVL('n'), JOBVR('n');
  CPPL_INT LDA(n), LDB(n), LDVL(1), LDVR(1), LWORK(8*n), INFO(1);
  double *BETA(new double[n]), *VL(NULL), *VR(NULL), *WORK(new double[LWORK]);
  dggev_(&JOBVL, &JOBVR, &n, array, &LDA, matB.array, &LDB, &wr[0], &wi[0], BETA, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] VL;
  delete [] VR;
  
  //// reforming ////
  for(CPPL_INT i=0; i<n; i++){
    wr[i]/=BETA[i];
    wi[i]/=BETA[i];
  }
  delete [] BETA; 
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate generalized eigenvalues and generalized right eigenvectors\n
  All of the arguments don't need to be initialized.
  wr, wi, vrr and vri are overwitten and become 
  real and imaginary part of generalized eigenvalue 
  and generalized right eigenvector, respectively. 
  This matrix and matB are also overwritten.
*/
inline CPPL_INT dgematrix::dggev(dgematrix& matB, 
                             std::vector<double>& wr, std::vector<double>& wi, 
                             std::vector<dcovector>& vrr, 
                             std::vector<dcovector>& vri)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n){
    ERROR_REPORT;
    std::cerr << "This matrix is not a square matrix." << std::endl
              << "This matrix is (" << m << "x" << n << ")." << std::endl;
    exit(1);
  }
  if(matB.m!=n || matB.n!=n){
    ERROR_REPORT;
    std::cerr << "The matrix B is not a square matrix having the same size as \"this\" matrix." << std::endl
              << "The B matrix is (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  wr.resize(n);
  wi.resize(n);
  vrr.resize(n);
  vri.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    vrr[i].resize(n);
    vri[i].resize(n);
  }
  
  dgematrix VR(n,n);
  char JOBVL('n'), JOBVR('V');
  CPPL_INT LDA(n), LDB(n), LDVL(1), LDVR(n), LWORK(8*n), INFO(1);
  double *BETA(new double[n]), *VL(NULL), *WORK(new double[LWORK]);
  dggev_(&JOBVL, &JOBVR, &n, array, &LDA, matB.array, &LDB, &wr[0], &wi[0], BETA, VL, &LDVL, VR.array, &LDVR, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] VL;
  
  //// reforming ////
  for(CPPL_INT i=0; i<n; i++){
    wr[i]/=BETA[i];
    wi[i]/=BETA[i];
  }
  delete [] BETA; 
  
  //// forming ////
  for(CPPL_INT j=0; j<n; j++){
    if(fabs(wi[j])<DBL_MIN){
      for(CPPL_INT i=0; i<n; i++){
        vrr[j](i) = VR(i,j);
        vri[j](i) = 0.0;
      }
    }
    else{
      for(CPPL_INT i=0; i<n; i++){
        vrr[j](i)   = VR(i,j);
        vri[j](i)   = VR(i,j+1);
        vrr[j+1](i) = VR(i,j);
        vri[j+1](i) =-VR(i,j+1);
      }
      j++;
    }
  }
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! calculate generalized eigenvalues and generalized left eigenvectors\n
  All of the arguments don't need to be initialized.
  wr, wi, vlr and vli are overwitten and become 
  real and imaginary part of generalized eigenvalue 
  and generalized left eigenvector, respectively. 
  This matrix and matB are also overwritten.
*/
inline CPPL_INT dgematrix::dggev(dgematrix& matB, 
                             std::vector<double>& wr, std::vector<double>& wi, 
                             std::vector<drovector>& vlr, 
                             std::vector<drovector>& vli)
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=n){
    ERROR_REPORT;
    std::cerr << "This matrix is not a square matrix." << std::endl
              << "This matrix is (" << m << "x" << n << ")." << std::endl;
    exit(1);
  }
  if(matB.m!=n || matB.n!=n){
    ERROR_REPORT;
    std::cerr << "The matrix B is not a square matrix having the same size as \"this\" matrix." << std::endl
              << "The B matrix is (" << matB.m << "x" << matB.n << ")." << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  wr.resize(n);
  wi.resize(n);
  vlr.resize(n);
  vli.resize(n);
  for(CPPL_INT i=0; i<n; i++){
    vlr[i].resize(n);
    vli[i].resize(n);
  }
  
  dgematrix VL(n,n);
  char JOBVL('V'), JOBVR('n');
  CPPL_INT LDA(n), LDB(n), LDVL(n), LDVR(1), LWORK(8*n), INFO(1);
  double *BETA(new double[n]), *VR(NULL), *WORK(new double[LWORK]);
  dggev_(&JOBVL, &JOBVR, &n, array, &LDA, matB.array, &LDB, &wr[0], &wi[0], BETA, VL.array, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
  delete [] WORK;
  delete [] VR;
  
  //// reforming ////
  for(CPPL_INT i=0; i<n; i++){
    wr[i]/=BETA[i];
    wi[i]/=BETA[i];
  }
  delete [] BETA; 
  
  //// forming ////
  for(CPPL_INT j=0; j<n; j++){
    if(fabs(wi[j])<DBL_MIN){
      for(CPPL_INT i=0; i<n; i++){
        vlr[j](i) = VL(i,j);
        vli[j](i) = 0.0;
      }
    }
    else{
      for(CPPL_INT i=0; i<n; i++){
        vlr[j](i)   = VL(i,j);
        vli[j](i)   =-VL(i,j+1);
        vlr[j+1](i) = VL(i,j);
        vli[j+1](i) = VL(i,j+1);
      }
      j++;
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
/*! compute the singular value decomposition (SVD)\n
  The argument is dgbmatrix S.
  S doesn't need to be initialized.
  S is overwitten and become singular values.
  This matrix also overwritten.
*/
inline CPPL_INT dgematrix::dgesvd
(
 dgbmatrix& S
 )
{CPPL_VERBOSE_REPORT;
  char JOBU ='N';
  char JOBVT ='N';
  // M
  // N
  // A
  CPPL_INT LDA =m;
  double* U =NULL;
  S.resize(m,n,0,0);
  CPPL_INT LDU =1;
  double* VT =NULL;
  CPPL_INT LDVT =1;
  CPPL_INT LWORK =std::max(3*std::min(m,n)+std::max(m,n),5*std::min(m,n));
  double *WORK =new double[LWORK];
  CPPL_INT INFO =1;
  
  dgesvd_(&JOBU, &JOBVT, &m, &n, array, &LDA, S.array, U, &LDU, VT, &LDVT, WORK, &LWORK, &INFO);
  delete [] WORK;
  
  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  }
  return INFO;
}

//=============================================================================
/*! compute the singular value decomposition (SVD)\n
  The arguments are dcocector S, dgematrix U and VT.
  All of them need not to be initialized.
  S, U and VT are overwitten and become singular values, 
  left singular vectors,
  and right singular vectors respectively.
  This matrix also overwritten.
*/
inline CPPL_INT dgematrix::dgesvd(dcovector& S, dgematrix& U, dgematrix& VT)
{CPPL_VERBOSE_REPORT;
  char JOBU('A'), JOBVT('A');
  CPPL_INT LDA(m), LDU(m), LDVT(n), LWORK(std::max(3*std::min(m,n)+std::max(m,n),5*std::min(m,n))), INFO(1);
  double *WORK(new double[LWORK]);
  S.resize(std::min(m,n));
  U.resize(LDU,m);
  VT.resize(LDVT,n);
  
  dgesvd_(&JOBU, &JOBVT, &m, &n, array, &LDA, S.array, U.array, &LDU, VT.array, &LDVT, WORK, &LWORK, &INFO);
  delete [] WORK;
  
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
/*! solve the linear equality-constrained least squares (LSE) problem\n
  Input matrix and vectors, B, c, and d, are overwitten.
  This matrix is also overwritten.
  The solution vector x is to be automatically resized.
*/
inline CPPL_INT dgematrix::dgglse
(
 dgematrix& B,
 dcovector& c,
 dcovector& d,
 dcovector& x
 )
{CPPL_VERBOSE_REPORT;
#ifdef  CPPL_DEBUG
  if(m!=c.l){
    ERROR_REPORT;
    std::cerr << "A.m and c.l should be the same." << std::endl
              << "Your input was A.m=" << m << " and c.l=" << c.l << std::endl;
    exit(1);
  }
  if(B.m!=d.l){
    ERROR_REPORT;
    std::cerr << "B.m and d.l should be the same." << std::endl
              << "Your input was B.m=" << B.m << " and d.l=" << d.l << std::endl;
    exit(1);
  }
  if( !(B.m<=n) || !(n<=m+B.m) ){
    ERROR_REPORT;
    std::cerr << "B.m<=A.n<=A.m+B.m should be satisfied." << std::endl
              << "Your input was B.m=" << B.m << ", A.n=" << n << ", and A.m+B.m=" << m+B.m << std::endl;
    exit(1);
  }
#endif//CPPL_DEBUG
  
  CPPL_INT lwork(-1), info(1);
  dcovector work(1);
  x.resize(n);
  
  //////// workspace query ////////
  CPPL_INT lda =std::max(CPPL_INT(1),m);
  CPPL_INT ldb =std::max(CPPL_INT(1),B.m);
  dgglse_(&m, &n, &B.m, array, &lda, B.array, &ldb, c.array, d.array, x.array, work.array, &lwork, &info);
  lwork =CPPL_INT(work(0));
  work.resize(lwork);
  info =1;
  
  //////// solve ////////
  dgglse_(&m, &n, &B.m, array, &lda, B.array, &ldb, c.array, d.array, x.array, work.array, &lwork, &info);
  work.clear();
  
  if(info!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. info = " << info << "." << std::endl;
  }
  return info;
}
