#include "gtest.h"
#include "fft.h"
#include <cstdio>
#include <cmath>
// #include <cassert>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);


double beta=10;
int N=1024;
double e0 = 0.4;
double tail = 0.77;
double acc_fft = 1e-3;
double acc_rounding = 1e-14;


void eval_G_tau(vector<double> &G_tau)
{
    G_tau.resize(N);
    for(int i=0; i<N; i++){
        double tau = beta * double(i) / double(N);
        G_tau[i] = -tail * exp(-e0*tau) / ( exp(-e0*beta) + 1.0 );
    }
}

void eval_G_iw(vector< complex<double> > &G_iw)
{
    G_iw.resize(N/2);
    for(int i=0; i<N/2; i++){
        double omega = double(2*i+1) * M_PI / beta;
        G_iw[i] = tail / (IMAG*omega - e0);
    }
}

template <typename T>
double diff_vector(const vector<T> &vec1, const vector<T> &vec2)
// double diff_vector(vector<T> vec1, vector<T> vec2)
{
    // assert(vec1.size() == vec2.size());
    double dif_max = 0;
    for(int i=0; i<vec1.size(); i++){
        dif_max = max( dif_max, abs(vec1[i] - vec2[i]) );
    }
    return dif_max;
}

// G(tau) -> G(iw)
TEST(FFT, tau2iw)
{
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact);

    vector< complex<double> > G_iw_fft;
    fft_fermion_tau2iw(G_tau_exact, G_iw_fft, beta, tail);
    double dif_iw = diff_vector(G_iw_exact, G_iw_fft);
    printf("dif_iw: %.3e\n", dif_iw);
    EXPECT_GT(acc_fft, dif_iw) << "Failed in G(tau) -> G(iw)";
}

// G(iw) -> G(tau)
TEST(FFT, iw2tau)
{
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact);

    vector<double> G_tau_fft;
    fft_fermion_iw2tau(G_tau_fft, G_iw_exact, beta, tail);
    double dif_tau = diff_vector(G_tau_exact, G_tau_fft);
    printf("dif_tau: %.3e\n", dif_tau);
    EXPECT_GT(acc_fft, dif_tau) << "Failed in G(iw) -> G(tau)";
}

// // G(iw) -> G(tau) -> G(iw)
TEST(FFT, iw2tau2iw)
{
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact);

    vector<double> G_tau_fft;
    vector< complex<double> > G_iw_fft;
    fft_fermion_iw2tau(G_tau_fft, G_iw_exact, beta, tail);
    fft_fermion_tau2iw(G_tau_fft, G_iw_fft, beta, tail);
    double dif_fft = diff_vector(G_iw_exact, G_iw_fft);
    printf("dif_fft: %.3e\n", dif_fft);
    EXPECT_GT(acc_rounding, dif_fft) << "Failed in G(iw) -> G(tau) -> G(iw)";
}
