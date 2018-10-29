#include "fft.h"
#include <cstdio>
#include <cmath>
#include <fftw3.h>
#include <complex>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);

static void fft_tau2iw(const std::vector<double> &g_tau, std::vector< std::complex<double> > &g_iw, const double beta)
{
    int N = g_tau.size();

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*N);
    p = fftw_plan_dft_1d(2*N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    for(int i=0; i<N; i++){
        in[i][0] = g_tau[i];
        in[i][1] = 0;
        in[i+N][0] = -g_tau[i];
        in[i+N][1] = 0;
    }

    fftw_execute(p);

    for(int i=0; i<N/2; i++){
        int n_odd = 2*i+1;
        g_iw[i] = complex<double>(out[n_odd][0], out[n_odd][1]) / double(2*N) * beta;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

static void fft_iw2tau(std::vector<double> &g_tau, const std::vector< std::complex<double> > &g_iw, const double beta)
{
    int N = g_iw.size()*2;

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*N);
    p = fftw_plan_dft_1d(2*N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for(int i=0; i<2*N; i++)  in[i][0] = in[i][1] = 0;

    for(int i=0; i<N/2; i++){
        int n_odd = 2*i+1;
        in[n_odd][0] = g_iw[i].real();
        in[n_odd][1] = g_iw[i].imag();
        in[2*N-n_odd][0] = g_iw[i].real();
        in[2*N-n_odd][1] = -g_iw[i].imag();
    }

    fftw_execute(p);

    for(int i=0; i<N; i++){
        g_tau[i] = out[i][0] / beta;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

void fft_fermion_tau2iw(const std::vector<double> &G_tau, std::vector< std::complex<double> > &G_iw, const double beta, const double tail)
{
    double N = G_tau.size();
    G_iw.resize(N/2);
    vector<double> G_tau_dif(N);

    for(int i=0; i<N; i++){
        G_tau_dif[i] = G_tau[i] + tail * 0.5;
    }

    fft_tau2iw(G_tau_dif, G_iw, beta);

    for(int i=0; i<N/2; i++){
        double iw = (double)(2*i+1) * M_PI / beta;
        G_iw[i] += tail / (IMAG*iw);
    }
}

void fft_fermion_iw2tau(std::vector<double> &G_tau, const std::vector< std::complex<double> > &G_iw, const double beta, const double tail)
{
    double N = G_iw.size()*2;
    G_tau.resize(N);
    vector< complex<double> > G_iw_dif(N/2);

    for(int i=0; i<N/2; i++){
        double iw = (double)(2*i+1) * M_PI / beta;
        G_iw_dif[i] = G_iw[i] - tail / (IMAG*iw);
    }

    fft_iw2tau(G_tau, G_iw_dif, beta);

    for(int i=0; i<N; i++){
        G_tau[i] -= tail * 0.5;
    }
}
