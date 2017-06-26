//=============================================================================
/*! solve A*x=b for real and unsymmetric matrix using Intel PARDISO.\n
  The argument is dcovector b.
  b is overwritten and become the solution x.
  A is not overwritten.
*/
inline CPPL_INT dgsmatrix::pardiso(dcovector& b) const
{CPPL_VERBOSE_REPORT;
#ifndef  __INTEL_COMPILER
  ERROR_REPORT;
  std::cerr << "dgsmatrix::pardiso is only for intel c++ compiler (icpc)." << std::endl;
  std::cerr << "Recompile your code with icpc to use pardiso." << std::endl;
  (void)b;
  exit(1);
#else  //__INTEL_COMPILER
  
  
  //#ifdef  CPPL_DEBUG
  if(m!=n || m!=b.l){
    ERROR_REPORT;
    std::cerr << "These matrix and vector cannot be solved." << std::endl
              << "Your input was (" << m << "x" << n << ") and (" << b.l << ")." << std::endl;
    exit(1);
  }
  //#endif//CPPL_DEBUG 
  
  //////// convert matrix storage into compressed sparse row format ////////
  //std::cerr << "converting" << std::endl;
  std::vector<double> a(data.size());
  std::vector<MKL_INT> ja(data.size());
  std::vector<MKL_INT> ia(m+1);
  ia[0] =0;
  MKL_INT k=0;
  for(CPPL_INT i=0; i<m; i++){
    //// make map ////
    const std::vector<CPPL_INT>::const_iterator rows_i_end =rows[i].end();
    std::map<CPPL_INT,CPPL_INT> jc;
    for(std::vector<CPPL_INT>::const_iterator rit=rows[i].begin(); rit!=rows_i_end; rit++){
      jc.insert( std::make_pair(data[*rit].j, *rit) );
    }
    //// assign ////
    const std::map<CPPL_INT,CPPL_INT>::const_iterator jc_end =jc.end();
    for(std::map<CPPL_INT,CPPL_INT>::const_iterator jcit=jc.begin(); jcit!=jc_end; jcit++){
      a[k] =data[(*jcit).second].v;
      ja[k] =MKL_INT((*jcit).first);
      k++;
    }
    ia[i+1] =k;
  }
  
  //////// pardisoinit ////////
  //std::cerr << "initializing" << std::endl;
  _MKL_DSS_HANDLE_t pt[64];
  MKL_INT mtype =11;//real unsymmetric
  MKL_INT iparm[64];
  PARDISOINIT(pt, &mtype, iparm);
  iparm[1] =3;//parallel fill-in reducing ordering
  iparm[23] =1;//use two-level scheduling factorization algorithm
  iparm[26] =0;//disable matrix checker
  iparm[34] =-1;//use zero-base array index
  
  //////// pardiso ////////
  //std::cerr << "solving" << std::endl;
  MKL_INT maxfct =1;
  MKL_INT mnum =1;
  MKL_INT phase =13;
  MKL_INT MKL_INT_n =MKL_INT(n);
  std::vector<MKL_INT> perm(n);
  MKL_INT nrhs =1;
  MKL_INT msglvl =0;
  dcovector x(b.l);
  MKL_INT error =1;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &MKL_INT_n, &a[0], &ia[0], &ja[0], &perm[0], &nrhs, iparm, &msglvl, b.array, x.array, &error);
  swap(b,x);//set b as x
  if(error!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. error = " << error << "." << std::endl;
  }
  
  //////// release memory ////////
  phase =-1;
  MKL_INT error2 =1;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &MKL_INT_n, &a[0], &ia[0], &ja[0], &perm[0], &nrhs, iparm, &msglvl, b.array, x.array, &error2);
  
  return error;
#endif //__INTEL_COMPILER
}
