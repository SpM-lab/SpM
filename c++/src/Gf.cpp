#include "Gf.h"
#include "fft.h"
#include <iostream>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);

void Gf::set_Gtau(const std::vector<double> &G_tau, const statistics stat, const double beta, const double tail)
{
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
    this->stat = stat;
    this->beta = beta;
    this->tail = tail;

    if(stat == FERMION){
        int n = G_iw.size();
        N = 2*n;
        giw.resize(N);
        for(int i=0; i<n; ++i){
            giw[i] = G_iw[i];
            giw[N-i-1] = conj(G_iw[i]);
        }
    }else{ // BOSON
        int n = G_iw.size();
        N = 2*n-1;
        giw.resize(N);
        giw[0] = G_iw[0];
        for(int i=1; i<n; ++i){
            giw[i] = G_iw[i];
            giw[N-i] = conj(G_iw[i]);
        }
    }
    gtau.clear();
    flag_pade = false;
}

void Gf::get_Gtau(std::vector<double> &G_tau)
{
    if( gtau.empty() )  compute_iw2tau();
    G_tau.assign(gtau.begin(), gtau.end());
}

void Gf::get_Giw(std::vector< std::complex<double> > &G_iw)
{
    if( giw.empty() )  compute_tau2iw();
    int n = (stat==FERMION ? N/2 : N/2+1 );
    G_iw.assign(giw.begin(), giw.begin()+n);
}

// void get_Gomega(std::vector< std::complex<double> > &G_omega, const double omega_min, const double omega_max, const int n)
// {
//
// }

std::complex<double> Gf::Gomega(const double omega)
{
    if( !flag_pade ){
        if( giw.empty() )  compute_tau2iw();

        int n = N/2;
        if(isodd(N)) ++n;
        vector< complex<double> > matsubara(n);
        vector< complex<double> > Giw(n);
        int offset = stat==BOSON?0:1;
        for(int i=0; i<n; i++){
            matsubara[i] = (2*i+offset) * M_PI / beta * IMAG;
            Giw[i] = giw[i];
        }
        pade.init(matsubara, Giw);
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
        fft_boson_tau2iw(gtau, giw, beta, tail);
    }
}
void Gf::compute_iw2tau()
{
    if( stat == FERMION ){
        fft_fermion_iw2tau(gtau, giw, beta, tail);
    }
    else{
        fft_boson_iw2tau(gtau, giw, beta, tail);
    }
}
