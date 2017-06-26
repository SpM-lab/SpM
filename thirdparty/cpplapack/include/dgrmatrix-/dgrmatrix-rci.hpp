//=============================================================================
/*! solve A*x=b for real and unsymmetric matrix using DFGMRES of Intel RCI ISS without preconditioning.\n
  The argument is dcovector b.
  b is overwritten and become the solution x.
  A is not overwritten.
*/
inline CPPL_INT dgrmatrix::dfgmres(dcovector& b, const double eps) const
{CPPL_VERBOSE_REPORT;
#ifndef  __INTEL_COMPILER
  ERROR_REPORT;
  std::cerr << "dgrmatrix::dfgmres is only for intel c++ compiler (icpc)." << std::endl;
  std::cerr << "Recompile your code with icpc to use dfgmres." << std::endl;
  (void)b;
  (void)eps;
  exit(1);
  
  
#else  //__INTEL_COMPILER is defined.
  //#ifdef  CPPL_DEBUG
  if(m!=n || m!=b.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << b.l << ")." << std::endl;
    exit(1);
  }
  //#endif//CPPL_DEBUG 
  
  //////////////// constant ////////////////
  MKL_INT itmax =1000;
  MKL_INT n_res =500;//number of iteration before restart
  MKL_INT N =MKL_INT(n);
  const dcovector b_orig =b;
  
  //////////////// dfgmresinit ////////////////
  dcovector x =dcovector(N).zero();//initial guess
  //dcovector x =b;//initial guess
  MKL_INT RCI_request;
  MKL_INT ipar[128];
  double dpar[128];
  std::vector<double> tmp((2*n_res+1)*N + (n_res*(n_res+9))/2 + 1);
  DFGMRES_INIT(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0]);
  if(RCI_request){
    ERROR_REPORT;
    std::cerr << "dfgmres_init failed. RCI_request=" << RCI_request << "." << std::endl;
    exit(1);
  }
  
  //////////////// dfgmres_check ////////////////
  //////// ipar ////////
  ipar[7] = 0;//if check iteration count
  ipar[8] = 0;//if enable residual stopping test
  ipar[9] = 1;//if enable user defined stopping test
  ipar[10] = 0;//if enable preconditioner
  ipar[11] = 1;//if enable check of the norm of the next generated vector automatically
  ipar[14] = n_res;//number of iteration before restart
  //////// dpar ////////
  ////dpar[0] =1e-99;//relative tolerance
  ////dpar[1] =0.;//absolute tolerance
  //////// call ////////
  DFGMRES_CHECK(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0]);
  if(RCI_request){
    ERROR_REPORT;
    std::cerr << "dfgmres_init failed. RCI_request=" << RCI_request << "." << std::endl;
    exit(1);
  }
  
  //////////////// dfgmres ////////////////
  char transa ='N';
  double* A =const_cast<double*>(&a[0]);
  MKL_INT* IA =const_cast<MKL_INT*>(&ia[0]);
  MKL_INT* JA =const_cast<MKL_INT*>(&ja[0]);
  MKL_INT itercount;
  dcovector x_tmp(N);
  size_t itc =0;
  while(1){
    //// call ////
    DFGMRES(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0]);
    
    //// RCI_request ////
    if(RCI_request==1){//continue
      MKL_DCSRGEMV(&transa, &N, A, IA, JA, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);//calc A*v
      itc++;
      //std::cerr << "A*v" << std::endl;
    }
    else if(RCI_request==2){//convergence check with count check
      if(itc%100==0){
        ipar[12] = 1;
        DFGMRES_GET(&N, x.array, x_tmp.array, &RCI_request, ipar, dpar, &tmp[0], &itercount);
        double norm_r =fabs(damax(b_orig-(*this)*x_tmp));
        std::cerr << "itc=" << itc << ": norm_r = " << norm_r << " (eps=" << eps << ")" << std::endl;
        if(norm_r<eps){
          std::cerr << "converged (norm_r<eps)" << std::endl;
          break;
        }
        if(itc>=itmax){
          std::cerr << "failed. (itc>=itmax)" << std::endl;
          break;
        }
      }
    }
    else if(RCI_request==0){//converged (ipar[11])
      std::cerr << "converged (RCI_request=0)" << std::endl;
      break;
    }
    else{//failed
      WARNING_REPORT;
      std::cerr << "dfgmres failed. RCI_request=" << RCI_request << "." << std::endl;
      break;
    }
  }
  //////// info ////////
  CPPL_INT info =RCI_request;
  
  //////////////// dfgmres_get ////////////////
  ipar[12] =0;
  DFGMRES_GET(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0], &itercount);
  
  //////////////// MKL_Free_Buffers ////////////////
  MKL_Free_Buffers();
  
  //////////////// swap x and b ////////////////
  swap(x, b);
  
  return info;
#endif //__INTEL_COMPILER
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
/*! solve A*x=b for real and unsymmetric matrix using DFGMRES with ILUT precondition of Intel RCI ISS.\n
  The argument is dcovector b.
  b is overwritten and become the solution x.
  A is not overwritten.
*/
inline CPPL_INT dgrmatrix::ilut_dfgmres(dcovector& b, const int maxfil, const double eps) const
{CPPL_VERBOSE_REPORT;
#ifndef  __INTEL_COMPILER
  ERROR_REPORT;
  std::cerr << "dgrmatrix::ilut_dfgmres is only for intel c++ compiler (icpc)." << std::endl;
  std::cerr << "Recompile your code with icpc to use ilut_dfgmres." << std::endl;
  (void)b;
  (void)maxfil;
  (void)eps;
  exit(1);
  
  
#else  //__INTEL_COMPILER is defined.
  //#ifdef  CPPL_DEBUG
  if(m!=n || m!=b.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << b.l << ")." << std::endl;
    exit(1);
  }
  //#endif//CPPL_DEBUG 
  
  //////////////// constants ////////////////
  MKL_INT itmax =500;
  MKL_INT n_res =500;//number of iteration before restart
  //MKL_INT n_res =150;//number of iteration before restart
  MKL_INT N =MKL_INT(n);
  MKL_INT MAXFIL =std::min(maxfil, N-1);//limitter for maxfil
  const dcovector b_orig =b;
  
  //////////////// dfgmresinit ////////////////
  dcovector x =dcovector(b.l).zero();//initial guess
  //dcovector x =b;//initial guess
  MKL_INT RCI_request;
  MKL_INT ipar[128];
  double dpar[128];
  std::vector<double> tmp((2*n_res+1)*N + (n_res*(n_res+9))/2 + 1);
  DFGMRES_INIT(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0]);
  if(RCI_request){
    ERROR_REPORT;
    std::cerr << "dfgmres_init failed. RCI_request=" << RCI_request << "." << std::endl;
    exit(1);
  }

  //////////////// dcsrilut ////////////////
  double* A =const_cast<double*>(&a[0]);
  MKL_INT* IA =const_cast<MKL_INT*>(&ia[0]);
  MKL_INT* JA =const_cast<MKL_INT*>(&ja[0]);
  std::vector<double> bilut((2*MAXFIL+1)*N - MAXFIL*(MAXFIL+1) + 1);
  std::vector<MKL_INT> ibilut(N+1);
  std::vector<MKL_INT> jbilut((2*MAXFIL+1)*N - MAXFIL*(MAXFIL+1) + 1);
  double tol =1e-3;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //double tol =1e-4;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //double tol =1e-6;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MKL_INT ierr;
  //ipar[30] =0;//dangerrous
  ipar[30] =1;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dpar[30] =tol;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //dpar[30] =1e-99;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DCSRILUT(&N, A, IA, JA, &bilut[0], &ibilut[0], &jbilut[0], &tol, &MAXFIL, ipar, dpar, &ierr);
  if(ierr){
    ERROR_REPORT;
    std::cerr << "dcsrilut failed. ierr=" << ierr << "." << std::endl;
    exit(1);
  }
  
  //////////////// dfgmres_check ////////////////
  //////// ipar ////////
  ipar[7] = 0;//if check iteration count
  ipar[8] = 0;//if enable residual stopping test
  ipar[9] = 1;//if enable user defined stopping test
  ipar[10] = 1;//if enable preconditioner
  ipar[11] = 1;//if enable check of the norm of the next generated vector automatically
  ipar[14] = n_res;//number of iteration before restart
  //////// dpar ////////
  ////dpar[0] =1e-1;//relative tolerance
  ////dpar[1] =eps;//absolute tolerance
  //////// call ////////
  DFGMRES_CHECK(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0]);
  if(RCI_request){
    ERROR_REPORT;
    std::cerr << "dfgmres_init failed. RCI_request=" << RCI_request << "." << std::endl;
    exit(1);
  }
  
  //////////////// dfgmres ////////////////
  char transa ='N';
  char uplo1   ='L';
  char transa1 ='N';
  char diag1   ='U';
  char uplo2   ='U';
  char transa2 ='N';
  char diag2   ='N';
  MKL_INT itercount;
  dcovector x_tmp(N);
  size_t itc =0;
  while(1){
    //// call ////
    DFGMRES(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0]);
    
    //// RCI_request ////
    if(RCI_request==1){//continue
      MKL_DCSRGEMV(&transa, &N, A, IA, JA, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);//calc A*v
      itc++;
    }
    else if(RCI_request==3){//precondition
      MKL_DCSRTRSV(&uplo1, &transa1, &diag1, &N, &bilut[0], &ibilut[0], &jbilut[0], &tmp[ipar[21]-1], x_tmp.array);
      MKL_DCSRTRSV(&uplo2, &transa2, &diag2, &N, &bilut[0], &ibilut[0], &jbilut[0], x_tmp.array, &tmp[ipar[22]-1]);
      //std::cout << "preconditioned" << std::endl;
    }
    else if(RCI_request==2){//convergence check with count check
      if(itc%100==0){
        ipar[12] = 1;
        DFGMRES_GET(&N, x.array, x_tmp.array, &RCI_request, ipar, dpar, &tmp[0], &itercount);
        double norm_r =fabs(damax(b_orig-(*this)*x_tmp));
        std::cerr << "itc=" << itc << ": norm_r = " << norm_r << " (eps=" << eps << ")" << std::endl;
        if(norm_r<eps){
          std::cerr << "converged (norm_r<eps)" << std::endl;
          break;
        }
        if(itc>=itmax){
          std::cerr << "failed. (itc>=itmax)" << std::endl;
          break;
        }
      }
    }
    else if(RCI_request==0){//converged (ipar[11])
      std::cerr << "converged (RCI_request=0)" << std::endl;
      break;
    }
    else{//failed
      WARNING_REPORT;
      std::cerr << "dfgmres failed. RCI_request=" << RCI_request << "." << std::endl;
      break;
    }
  }
  //////// info ////////
  CPPL_INT info =RCI_request;
  
  //////////////// dfgmres_get ////////////////
  ipar[12] =0;
  DFGMRES_GET(&N, x.array, b.array, &RCI_request, ipar, dpar, &tmp[0], &itercount);
  
  //////////////// MKL_Free_Buffers ////////////////
  MKL_Free_Buffers();
  
  //////////////// swap x and b ////////////////
  swap(x, b);
  
  return info;
#endif //__INTEL_COMPILER
}
