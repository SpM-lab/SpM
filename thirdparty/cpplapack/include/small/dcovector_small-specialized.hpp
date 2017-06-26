//=============================================================================
/*! calculate vector product only for 2D vector */
inline double operator/(const dcovec2& A, const dcovec2& B)
{CPPL_VERBOSE_REPORT;
  return A(0)*B(1) -A(1)*B(0);
}

//=============================================================================
/*! convert 2D vector to theta */
inline double v2t(const dcovec2& v)
{CPPL_VERBOSE_REPORT;
  return std::atan2(v(1),v(0));
}

//=============================================================================
/*! rotate 2D vector by theta [rad] */
inline dcovec2 rotate(const dcovec2& v, const double& t)
{CPPL_VERBOSE_REPORT;
  dcovec2 w;
  w(0) =v(0)*std::cos(t) -v(1)*std::sin(t);
  w(1) =v(0)*std::sin(t) +v(1)*std::cos(t);
  return w;
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
/*! calculate vector product only for 3D vector */
inline dcovec3 operator/(const dcovec3& A, const dcovec3& B)
{CPPL_VERBOSE_REPORT;
  dcovec3 C;
  C(0) =A(1)*B(2) -A(2)*B(1);
  C(1) =A(2)*B(0) -A(0)*B(2);
  C(2) =A(0)*B(1) -A(1)*B(0);
  return C;
}

//=============================================================================
/*! calculate vector product only for 3D vector */
inline dcovec3 operator/=(dcovec3& A, const dcovec3& B)
{CPPL_VERBOSE_REPORT;
  A =A/B;
  return A;
}

//=============================================================================
/*! make quaternion from imag vector and real value */
inline dquater vr2q(const dcovec3& v, const double& r)
{CPPL_VERBOSE_REPORT;
  return dquater(v(0),v(1),v(2),r);
}

//=============================================================================
/*! make quaternion from directional vector and rotational angle */
inline dquater vt2q(const dcovec3& v, const double& theta)
{CPPL_VERBOSE_REPORT;
  return vr2q( v/(nrm2(v)+DBL_MIN)*std::sin(0.5*theta), std::cos(0.5*theta) );
}

//=============================================================================
/*! rotate 3D vector by quaternion */
inline dcovec3 rotate(const dcovec3& v, const dquater& q)
{CPPL_VERBOSE_REPORT;
  return imag( q*vr2q(v,0.)*conj(q) );
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
/*! conjuction */
inline dquater conj(const dquater& q)
{CPPL_VERBOSE_REPORT;
  return dquater(-q(0),-q(1),-q(2), q(3));
}

//=============================================================================
/*! imag */
inline dcovec3 imag(const dquater& q)
{CPPL_VERBOSE_REPORT;
  return dcovec3(q(0),q(1),q(2));
}

//=============================================================================
/*! inverse */
inline dquater inv(const dquater& q)
{CPPL_VERBOSE_REPORT;
  return conj(q)/std::pow(nrm2(q),2);
}

//=============================================================================
/*! dquater*dquater operator */
inline dquater operator*(const dquater& q1, const dquater& q2)
{CPPL_VERBOSE_REPORT;
  return dquater(q1(3)*q2(0) +q1(0)*q2(3) +q1(1)*q2(2) -q1(2)*q2(1),
                 q1(3)*q2(1) -q1(0)*q2(2) +q1(1)*q2(3) +q1(2)*q2(0),
                 q1(3)*q2(2) +q1(0)*q2(1) -q1(1)*q2(0) +q1(2)*q2(3),
                 q1(3)*q2(3) -q1(0)*q2(0) -q1(1)*q2(1) -q1(2)*q2(2) );
}

//=============================================================================
/*! dquater/dquater operator */
inline dquater operator/(const dquater& q1, const dquater& q2)
{CPPL_VERBOSE_REPORT;
  return q1*inv(q2);
}

//=============================================================================
/*! dquater*=dquater operator */
inline dquater operator*=(dquater& q1, const dquater& q2)
{CPPL_VERBOSE_REPORT;
  q1 =q1*q2;
  return q1;
}

//=============================================================================
/*! dquater/=dquater operator */
inline dquater operator/=(dquater& q1, const dquater& q2)
{CPPL_VERBOSE_REPORT;
  q1 =q1/q2;
  return q1;
}

//=============================================================================
/*! return vector from quaternion (|vector|=theta) */
inline dcovec3 q2vt(const dquater& q)
{CPPL_VERBOSE_REPORT;
  double sin_theta_half;
  double theta( 2.*std::acos(q(3)) );
  
  if(theta<M_PI){
    sin_theta_half =std::sin(0.5*theta);
  }
  else{
    theta -=2.*M_PI;
    sin_theta_half =-std::sin(0.5*theta);
  }
  
  return dcovec3( theta*q(0)/sin_theta_half,
                  theta*q(1)/sin_theta_half,
                  theta*q(2)/sin_theta_half );
}

//=============================================================================
/*! return rotational matrix made of quaternion */
inline dgemat3 q2m(const dquater& q)
{CPPL_VERBOSE_REPORT;
  dquater cq( conj(q) );
  dquater X( dquater(+q(3),+q(2),-q(1),-q(0))*cq );
  dquater Y( dquater(-q(2),+q(3),+q(0),-q(1))*cq );
  dquater Z( dquater(+q(1),-q(0),+q(3),-q(2))*cq );
  dgemat3 mat;
  mat(0,0)=X(0); mat(0,1)=Y(0); mat(0,2)=Z(0);
  mat(1,0)=X(1); mat(1,1)=Y(1); mat(1,2)=Z(1);
  mat(2,0)=X(2); mat(2,1)=Y(2); mat(2,2)=Z(2);
  return mat;
}
