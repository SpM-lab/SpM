/*!
 \file pade_c.cpp
 \brief Definition of class Pade
*/

#include "pade.h"

using namespace std;


Pade::Pade(vector<complex<double> > &z, vector<complex<double> > &u)
{
	init(z, u);
}

void Pade::init(vector<complex<double> > &z, vector<complex<double> > &u)
{
	m_n = min(z.size(), u.size());
	m_z = z;  // copy
	m_a.resize(m_n);

	// g_n(z_i)
	vector<complex<double> > *g0 = new vector<complex<double> >(u);  // copy
	vector<complex<double> > *g1 = new vector<complex<double> >(m_n);

	m_a[0] = u[0];
	for(int n=1; n<m_n; n++){
		// g0[i] = g_{n-1}[i]
		// g1[i] = g_{n}[i]
		for(int i=n; i<m_n; i++)  (*g1)[i] = ( (*g0)[n-1] / (*g0)[i] - 1.0 ) / (z[i] - z[n-1]);
		m_a[n] = (*g1)[n];  // g[n][n]
		swap(g0, g1);
	}

	delete g0, g1;

	// truncate the coefficients m_a if m_a[n]=0
	for(int i=0; i<m_n; i++){
		// printf("%d %.15e %.15e\n", i, m_a[i].real(), m_a[i].imag());
		if( abs(m_a[i]) < 1e-12 ){  // if zero
			// printf("m_n=%d\n", m_n);
			m_n = i;
			break;
		}
	}
}

complex<double> Pade::eval(complex<double> w)
{
	complex<double> A0 = 0.0, A1 = m_a[0], A2, B0 = 1.0, B1 = 1.0, B2;

	for(int i=1; i<m_n; i++){  // (N-1)/2
		A2 = A1 + ( w - m_z[i-1] ) * m_a[i] * A0;
		B2 = B1 + ( w - m_z[i-1] ) * m_a[i] * B0;

		A1 /= B2;
		A2 /= B2;
		B1 /= B2;
		// B2 /= B2;
		B2 = 1.;

		A0 = A1;
		A1 = A2;
		B0 = B1;
		B1 = B2;
	}

	return( A2 / B2 );
}
