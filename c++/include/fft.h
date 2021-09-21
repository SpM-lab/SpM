#ifndef _FFT_H
#define _FFT_H

#include <vector>
#include <complex>

inline bool isodd(int n) { return (n/2)*2 != n; }
inline bool iseven(int n){ return !isodd(n); }

void fft_fermion_tau2iw(const std::vector<double> &G_tau, std::vector< std::complex<double> > &G_iw, const double beta, const double tail=1.);

void fft_fermion_iw2tau(std::vector<double> &G_tau, const std::vector< std::complex<double> > &G_iw, const double beta, const double tail=1.);

void fft_boson_tau2iw(const std::vector<double> &G_tau, std::vector< std::complex<double> > &G_iw, const double beta, const double tail=0.);

void fft_boson_iw2tau(std::vector<double> &G_tau, const std::vector< std::complex<double> > &G_iw, const double beta, const double tail=0.);

#endif // _FFT_H
