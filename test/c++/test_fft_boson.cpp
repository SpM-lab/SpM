#include "gtest.h"
#include "fft.h"
#include <cstdio>
#include <cmath>
#include <iostream>
// #include <cassert>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);


double beta = 10.0;
double e0 = 0.4;
double tail = 0.77;
// double tail = 1.0;
double acc_fft = 1e-3;
double acc_rounding = 1e-14;

void eval_G_tau(vector<double> &G_tau, int N)
{
    G_tau.resize(N);
    for(int i=0; i<N; i++){
        double tau = i*beta/N;
        G_tau[i] = tail * exp(-e0*tau) / (exp(-e0*beta)-1.0);
    }
}

// void eval_G_iw(vector< complex<double> > &G_iw, int N)
// {
//     const double bN = beta/N;
//     G_iw.resize(N);
//     G_iw[0] = -bN*tail/(1.0-exp(-beta*e0/N));
//     for(int i=1; i<N/2; i++){
//         double omega = 2*i* M_PI / beta;
//         G_iw[i] = -bN*tail/(1.0-exp(-bN*(e0-IMAG*omega)));
//         G_iw[N-i] = conj(G_iw[i]);
//     }
//     if(iseven(N)){
//         double omega = N * M_PI / beta;
//         G_iw[N/2] = -bN*tail/(1.0-exp(-bN*(e0-IMAG*omega)));
//     }else{
//         double omega = (N-1) * M_PI / beta;
//         G_iw[N/2] = -bN*tail/(1.0-exp(-bN*(e0-IMAG*omega)));
//         G_iw[N/2+1] = conj(G_iw[N/2]);
//     }
// }
void eval_G_iw(vector< complex<double> > &G_iw, int N)
{
    G_iw.resize(N);
    G_iw[0] = tail / (-e0);
    for(int i=1; i<N/2; i++){
        double omega = 2*i* M_PI / beta;
        G_iw[i] = tail / (IMAG*omega - e0);
        G_iw[N-i] = conj(G_iw[i]);
    }
    if(iseven(N)){
        double omega = N * M_PI / beta;
        G_iw[N/2] = tail / (IMAG*omega - e0);
    }else{
        double omega = (N-1) * M_PI / beta;
        G_iw[N/2] = tail / (IMAG*omega - e0);
        G_iw[N/2+1] = conj(G_iw[N/2]);
    }
}

template <typename T>
double diff_vector(const vector<T> &vec1, const vector<T> &vec2)
// double diff_vector(vector<T> vec1, vector<T> vec2)
{
    assert(vec1.size() == vec2.size());
    double dif_max = 0;
    for(int i=0; i<vec1.size(); i++){
        dif_max = max( dif_max, abs(vec1[i] - vec2[i]) );
        // std::cout << i << " " << vec1[i] << " " << vec2[i] << " " << vec1[i]-vec2[i] << " " << dif_max << std::endl;
    }
    return dif_max;
}

// G(tau) -> G(iw)
TEST(FFT, tau2iw_even_N)
{
    const int N = 1024;
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact, N);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact, N);

    vector< complex<double> > G_iw_fft;
    fft_boson_tau2iw(G_tau_exact, G_iw_fft, beta, tail);
    double dif_iw = diff_vector(G_iw_exact, G_iw_fft);
    printf("dif_iw: %.3e\n", dif_iw);
    EXPECT_GT(acc_fft, dif_iw) << "Failed in G(tau) -> G(iw)";
}

// G(tau) -> G(iw)
TEST(FFT, tau2iw_odd_N)
{
    const int N = 1025;
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact, N);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact, N);

    vector< complex<double> > G_iw_fft;
    fft_boson_tau2iw(G_tau_exact, G_iw_fft, beta, tail);
    double dif_iw = diff_vector(G_iw_exact, G_iw_fft);
    printf("dif_iw: %.3e\n", dif_iw);
    EXPECT_GT(acc_fft, dif_iw) << "Failed in G(tau) -> G(iw)";
}

// G(iw) -> G(tau)
TEST(FFT, iw2tau_even_N)
{
    const int N = 1024;
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact, N);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact, N);

    vector<double> G_tau_fft;
    fft_boson_iw2tau(G_tau_fft, G_iw_exact, beta, tail);
    double dif_tau = diff_vector(G_tau_exact, G_tau_fft);
    printf("dif_tau: %.3e\n", dif_tau);
    EXPECT_GT(acc_fft, dif_tau) << "Failed in G(iw) -> G(tau)";
}

// G(iw) -> G(tau)
TEST(FFT, iw2tau_odd_N)
{
    const int N = 1025;
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact, N);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact, N);

    vector<double> G_tau_fft;
    fft_boson_iw2tau(G_tau_fft, G_iw_exact, beta, tail);
    double dif_tau = diff_vector(G_tau_exact, G_tau_fft);
    printf("dif_tau: %.3e\n", dif_tau);
    EXPECT_GT(acc_fft, dif_tau) << "Failed in G(iw) -> G(tau)";
}

/*
 * In the bosonic case, some information is dropped for even N.
 * Fortunately, in our usage in Gf, this problem does not appear,
 * and so disable these tests.
 */


// G(tau) -> G(iw) -> G(tau)
// TEST(FFT, tau2iw2tau_even_N)
// {
//     const int N = 1024;
//     vector<double> G_tau_exact;
//     eval_G_tau(G_tau_exact, N);
//
//     vector< complex<double> > G_iw_exact;
//     eval_G_iw(G_iw_exact, N);
//
//     vector<double> G_tau_fft;
//     vector< complex<double> > G_iw_fft;
//     fft_boson_tau2iw(G_tau_exact, G_iw_fft, beta, tail);
//     fft_boson_iw2tau(G_tau_fft, G_iw_fft, beta, tail);
//     double dif_fft = diff_vector(G_tau_exact, G_tau_fft);
//     printf("dif_fft: %.3e\n", dif_fft);
//     EXPECT_GT(acc_rounding, dif_fft) << "Failed in G(tau) -> G(iw) -> G(tau)";
// }

// G(iw) -> G(tau) -> G(iw)
// TEST(FFT, iw2tau2iw_even_N)
// {
//     const int N = 1024;
//     vector<double> G_tau_exact;
//     eval_G_tau(G_tau_exact, N);
//
//     vector< complex<double> > G_iw_exact;
//     eval_G_iw(G_iw_exact, N);
//
//     vector<double> G_tau_fft;
//     vector< complex<double> > G_iw_fft;
//     fft_boson_iw2tau(G_tau_fft, G_iw_exact, beta, tail);
//     fft_boson_tau2iw(G_tau_fft, G_iw_fft, beta, tail);
//     double dif_fft = diff_vector(G_iw_exact, G_iw_fft);
//     printf("dif_fft: %.3e\n", dif_fft);
//     EXPECT_GT(acc_rounding, dif_fft) << "Failed in G(iw) -> G(tau) -> G(iw)";
// }

// G(tau) -> G(iw) -> G(tau)
TEST(FFT, tau2iw2tau_odd_N)
{
    const int N = 1025;
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact, N);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact, N);

    vector<double> G_tau_fft;
    vector< complex<double> > G_iw_fft;
    fft_boson_tau2iw(G_tau_exact, G_iw_fft, beta, tail);
    fft_boson_iw2tau(G_tau_fft, G_iw_fft, beta, tail);
    double dif_fft = diff_vector(G_tau_exact, G_tau_fft);
    printf("dif_fft: %.3e\n", dif_fft);
    EXPECT_GT(acc_rounding, dif_fft) << "Failed in G(tau) -> G(iw) -> G(tau)";
}

// G(iw) -> G(tau) -> G(iw)
TEST(FFT, iw2tau2iw_odd_N)
{
    const int N = 1025;
    vector<double> G_tau_exact;
    eval_G_tau(G_tau_exact, N);

    vector< complex<double> > G_iw_exact;
    eval_G_iw(G_iw_exact, N);

    vector<double> G_tau_fft;
    vector< complex<double> > G_iw_fft;
    fft_boson_iw2tau(G_tau_fft, G_iw_exact, beta, tail);
    fft_boson_tau2iw(G_tau_fft, G_iw_fft, beta, tail);
    double dif_fft = diff_vector(G_iw_exact, G_iw_fft);
    printf("dif_fft: %.3e\n", dif_fft);
    EXPECT_GT(acc_rounding, dif_fft) << "Failed in G(iw) -> G(tau) -> G(iw)";
}
