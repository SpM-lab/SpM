#include "fft.h"
#include <cstdio>
#include <cmath>
#include <fftw3.h>
#include <complex>
#include <iostream>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);

void fft_fermion_tau2iw(const std::vector<double> &G_tau,
                        std::vector< std::complex<double> > &G_iw,
                        const double beta,
                        const double tail)
{
    // G_tau excludes G(beta)

    const int N = G_tau.size();
    const int N2 = 2*N;
    G_iw.resize(N);

    fftw_complex* in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N2));
    fftw_complex* out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N2));
    fftw_plan p = fftw_plan_dft_1d(N2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // anti-periodic
    for(int i=0; i<N; i++){
        in[i][0] = G_tau[i] + tail*0.5;
        in[i][1] = 0.0;
        in[i+N][0] = -in[i][0];
        in[i+N][1] = 0.0;
    }

    fftw_execute(p);

    const double coeff = beta/N2;
    const double pib = M_PI/beta;
    for(int i=0; i<N/2; i++){
        const int n = 2*i+1;
        double iw = n*pib;
        G_iw[i] = coeff*std::complex<double>(out[n][0], out[n][1]);
        G_iw[i] += tail / (IMAG*iw);
        G_iw[N-i-1] = conj(G_iw[i]);
    }
    if(isodd(N)){
        G_iw[N/2] = coeff*std::complex<double>(out[N][0], out[N][1]);
        G_iw[N/2] += tail / (pib*N*IMAG);
    }

    fftw_free(out);
    fftw_free(in);
    fftw_destroy_plan(p);
}

void fft_fermion_iw2tau(std::vector<double> &G_tau,
                        const std::vector< std::complex<double> > &G_iw,
                        const double beta,
                        const double tail)
{
    int N = G_iw.size();
    const int N2 = 2*N;
    G_tau.resize(N);

    fftw_complex* in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N2));
    fftw_complex* out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N2));
    fftw_plan p = fftw_plan_dft_1d(N2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


    for(int i=0; i<N2; ++i){
        in[i][0] = in[i][1] = 0.0;
    }

    for(int i=0; i<N/2; ++i){
        const int n = 2*i+1;
        const double iw = n * M_PI / beta;
        in[n][0] = G_iw[i].real();
        in[n][1] = G_iw[i].imag()+tail/iw;
        in[N2-n][0] = in[n][0];
        in[N2-n][1] = -in[n][1];
    }
    if(isodd(N)){
      in[N][0] = G_iw[N/2].real();
      in[N][1] = G_iw[N/2].real()+beta*tail/(M_PI*N);
    }

    fftw_execute(p);

    for(int i=0; i<N; i++){
        G_tau[i] = out[i][0]/beta - tail*0.5;
    }

    fftw_destroy_plan(p);
    fftw_free(out);
    fftw_free(in);
}

static void gen_g0_iw_boson(std::vector< std::complex<double> > &g0_iw, const double beta, const double tail)
{
    int N = g0_iw.size();
    {
        // g0_iw[0] = tail * beta / 2.;
        g0_iw[0] = 0;
        for(int i=1; i<N/2; i++){
            double omega = 2*i* M_PI / beta;
            g0_iw[i] = tail / (IMAG*omega);
            g0_iw[N-i] = conj(g0_iw[i]);
        }
        if(iseven(N)){
            double omega = N * M_PI / beta;
            g0_iw[N/2] = tail / (IMAG*omega);
        }else{
            double omega = (N-1) * M_PI / beta;
            g0_iw[N/2] = tail / (IMAG*omega);
            g0_iw[N/2+1] = conj(g0_iw[N/2]);
        }
    }
}

void fft_boson_tau2iw(const std::vector<double> &G_tau,
                      std::vector< std::complex<double> > &G_iw,
                      const double beta,
                      const double tail)
{
    // G_tau excludes G(beta)
    int N = G_tau.size();
    G_iw.resize(N);

    fftw_complex *in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex *out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));
    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    for(int i=0; i<N; i++){
        // TODO: check sign
        // in[i][0] = -G_tau[i];
        in[i][0] = G_tau[i] - tail * ( double(i) / double(N) - 0.5 );  // tail*(tau/beta-1/2)
        in[i][1] = 0.0;
    }

    fftw_execute(p);

    std::vector< complex<double> > g0_iw(N);  // tail
    gen_g0_iw_boson(g0_iw, beta, tail);

    for(int i=0; i<N; i++){
        // const double omega = 2*M_PI*i/beta;
        G_iw[i] = complex<double>(out[i][0], out[i][1]) * (beta/N) + g0_iw[i];
    }

    fftw_destroy_plan(p);
    fftw_free(out);
    fftw_free(in);
}

void fft_boson_iw2tau(std::vector<double> &G_tau, const std::vector< std::complex<double> > &G_iw, const double beta, const double tail)
{
    int N = G_iw.size();
    G_tau.resize(N);

    fftw_complex *in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));
    fftw_complex *out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N));
    fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    std::vector< complex<double> > g0_iw(N);  // tail
    gen_g0_iw_boson(g0_iw, beta, tail);

    for(int i=0; i<N; i++){
        in[i][0] = G_iw[i].real() - g0_iw[i].real();
        in[i][1] = G_iw[i].imag() - g0_iw[i].imag();
    }

    fftw_execute(p);

    for(int i=0; i<N; i++){
        // TODO: check sign
        // G_tau[i] = -out[i][0] / beta;
        G_tau[i] = out[i][0] / beta + tail * ( double(i) / double(N) - 0.5 );  // tail*(tau/beta-1/2)
    }

    fftw_destroy_plan(p);
    fftw_free(out);
    fftw_free(in);
}
