
#ifndef _GF_H
#define _GF_H

#include <vector>
#include <complex>
#include "pade.h"

class Gf{
public:
    Gf(): flag_pade(false) {};

    enum statistics{ FERMION, BOSON };

    // initialize by G(tau). tau=[0:beta), N points (N=2^integer is recommended)
    void set_Gtau(const std::vector<double> &G_tau, const statistics stat, const double beta, const double tail=1.);

    // initialize by G(iw_n). n=[0:N/2)
    void set_Giw(const std::vector< std::complex<double> > &G_iw, const statistics stat, const double beta, const double tail=1.);

    int get_size() { return N; };
    double get_tail() { return tail; };
    double get_beta() { return beta; };
    statistics get_statistics() { return stat; };

    // get G(tau)
    void get_Gtau(std::vector<double> &G_tau);

    // get G(iw)
    void get_Giw(std::vector< std::complex<double> > &G_iw);

    // get G(w). w=[omega_min:omega_max], n-points
    // void get_Gomega(std::vector< std::complex<double> > &G_omega, const double omega_min, const double omega_max, const int n);

    // return G(omega)
    std::complex<double> Gomega(const double omega);

    // return rho(omega)
    double rho(const double omega);

private:
    int N;
    double beta;
    double tail;
    statistics stat;
    std::vector<double> gtau;  // size N
    std::vector< std::complex<double> > giw;  // size N/2

    bool flag_pade;
    Pade pade;

    void compute_tau2iw();  // G(tau) --> G(iw)
    void compute_iw2tau();  // G(iw) --> G(tau)
};

#endif // _GF_H
