//=============================================================================
/*! calculate determinant */
inline double det(const dsymat2& A)
{CPPL_VERBOSE_REPORT;
  return A(0,0)*A(1,1) -A(1,0)*A(1,0);
}

//=============================================================================
/*! calculate inverse */
inline dsymat2 inv(const dsymat2& A)
{CPPL_VERBOSE_REPORT;
  const double Adet( det(A) );
  dsymat2 Ainv;
  Ainv(0,0)= A(1,1)/Adet;
  Ainv(1,0)=-A(1,0)/Adet;  Ainv(1,1)= A(0,0)/Adet;
  return Ainv;
}

//=============================================================================
/*! rotate 2D matrix by rotational angle */
inline dsymat2 rotate(const dsymat2& m, const double& theta)
{CPPL_VERBOSE_REPORT;
  double c(std::cos(theta)), s(std::sin(theta));
  double cc(c*c), cs(c*s), ss(s*s);
  dsymat2 mat;
  mat(0,0) =m(0,0)*cc -2.*m(1,0)*cs       +m(1,1)*ss;
  mat(1,0) =m(1,0)*cc +(m(0,0)-m(1,1))*cs -m(1,0)*ss;
  mat(1,1) =m(1,1)*cc +2.*m(1,0)*cs       +m(0,0)*ss;
  return mat;
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
/*! calculate determinant */
inline double det(const dsymat3& A)
{CPPL_VERBOSE_REPORT;
  return
    +A(0,0)*A(1,1)*A(2,2) -A(0,0)*A(2,1)*A(2,1)
    +A(1,0)*A(2,1)*A(2,0) -A(1,0)*A(1,0)*A(2,2)
    +A(2,0)*A(1,0)*A(2,1) -A(2,0)*A(1,1)*A(2,0);
}

//=============================================================================
/*! calculate inverse */
inline dsymat3 inv(const dsymat3& A)
{CPPL_VERBOSE_REPORT;
  const double Adet( det(A) );
  dsymat3 Ainv;
  Ainv(0,0) =(A(1,1)*A(2,2)-A(2,1)*A(2,1))/Adet;
  Ainv(1,0) =(A(2,1)*A(2,0)-A(1,0)*A(2,2))/Adet;
  Ainv(1,1) =(A(0,0)*A(2,2)-A(2,0)*A(2,0))/Adet;
  Ainv(2,0) =(A(1,0)*A(2,1)-A(1,1)*A(2,0))/Adet;
  Ainv(2,1) =(A(1,0)*A(2,0)-A(0,0)*A(2,1))/Adet;
  Ainv(2,2) =(A(0,0)*A(1,1)-A(1,0)*A(1,0))/Adet;
  return Ainv;
}

//=============================================================================
/*! rotate 3D matrix by quaternion */
inline dsymat3 rotate(const dsymat3& m, const dquater& q)
{CPPL_VERBOSE_REPORT;
  dgemat3 R =q2m(q);
  dgemat3 Rm =R*m;
  
  dsymat3 RmRT;//not dgemat3
  RmRT.zero();
  for(CPPL_INT i=0; i<3; i++){
    for(CPPL_INT j=0; j<=i; j++){
      for(CPPL_INT k=0; k<3; k++){
        RmRT(i,j) +=Rm(i,k)*R(j,k);
      }
    }
  }
  return RmRT;
}
