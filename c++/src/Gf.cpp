#include "Gf.h"
#include "fft.h"
#include <iostream>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);

void Gf::set_Gtau(const std::vector<double> &G_tau, const statistics stat, const double beta, const double tail)
{
    if( stat == BOSON ){
        // TODO
        cout << "BOSON not yet implemented" << endl;
        exit(0);
    }
    gtau = G_tau;  // copy
    N = G_tau.size();
    this->stat = stat;
    this->beta = beta;
    this->tail = tail;

    giw.clear();
    flag_pade = false;
}

// initialize by G(iw_n). n=[0:N/2)
void Gf::set_Giw(const std::vector< std::complex<double> > &G_iw, const statistics stat, const double beta, const double tail)
{
    if( stat == BOSON ){
        // TODO
        cout << "BOSON not yet implemented" << endl;
        exit(0);
    }
    giw = G_iw;  // copy
    N = G_iw.size()*2;
    this->stat = stat;
    this->beta = beta;
    this->tail = tail;

    gtau.clear();
    flag_pade = false;
}

void Gf::get_Gtau(std::vector<double> &G_tau)
{
    if( gtau.empty() )  compute_iw2tau();
    G_tau = gtau;  // copy
}

void Gf::get_Giw(std::vector< std::complex<double> > &G_iw)
{
    if( giw.empty() )  compute_tau2iw();
    G_iw = giw;  // copy
}

// void get_Gomega(std::vector< std::complex<double> > &G_omega, const double omega_min, const double omega_max, const int n)
// {
//
// }

std::complex<double> Gf::Gomega(const double omega)
{
    if( !flag_pade ){
        if( giw.empty() )  compute_tau2iw();

        vector< complex<double> > matsubara(N);
        // TODO: boson
        for(int i=0; i<N; i++)  matsubara[i] = double(2*i+1) * M_PI / beta * IMAG;
        pade.init(matsubara, giw);
        flag_pade = true;
    }
    return pade.eval(omega);
}

double Gf::rho(const double omega)
{
    return -imag(Gomega(omega))/M_PI;
}

void Gf::compute_tau2iw()
{
    if( stat == FERMION ){
        fft_fermion_tau2iw(gtau, giw, beta, tail);
    }
    else{
        // TODO: BOSON
    }
}
void Gf::compute_iw2tau()
{
    if( stat == FERMION ){
        fft_fermion_iw2tau(gtau, giw, beta, tail);
    }
    else{
        // TODO: BOSON
    }
}
