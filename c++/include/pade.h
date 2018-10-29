/*!
 \file pade_c.h
 \brief Declaration of class Pade
*/

#include <complex>
#include <vector>

//! Analytic continuation by the Pade approximation

//! HOW TO USE:
//!  Pade obj(z, u);  ||  Pade obj;  obj.init(z,u);
//!  for(w){
//!    g = obj.eval(w);
//!  }
class Pade{
public:
	Pade() {};
	Pade(std::vector<std::complex<double> > &z, std::vector<std::complex<double> > &u);
	//! initialize a rational function u_{pade}(w) using N values of (u_i, z_i)
	void init(std::vector<std::complex<double> > &z, std::vector<std::complex<double> > &u);
	//! evaluate u_{pade}(w)
	std::complex<double> eval(std::complex<double> w);
private:
	int m_n;  // size
	//! coefficient set of the approximated function 'u_{pade}(z)'
	std::vector<std::complex<double> > m_a;
	std::vector<std::complex<double> > m_z;
};
