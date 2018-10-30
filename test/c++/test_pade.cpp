#include "gtest.h"
#include "Gf.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iterator>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);

template <typename T>
void read_file_to_vector(const char *filename, vector<T> &out)
{
    std::ifstream fin(filename);
    if( fin.fail() ){
        std::cerr << "file '" << filename << "' not exists" << std::endl;
        exit(1);
    }
    std::copy(
        std::istream_iterator<T>(fin),
        std::istream_iterator<T>(),
        std::back_inserter(out));
    fin.close();
}

template <typename T>
double diff_vector(const vector<T> &vec1, const vector<T> &vec2)
{
    // assert(vec1.size() == vec2.size());
    double dif_max = 0;
    for(int i=0; i<vec1.size(); i++){
        dif_max = max( dif_max, abs(vec1[i] - vec2[i]) );
    }
    return dif_max;
}



int N=256;
double beta=15;

double wmax = 1.5;
double wmin = -wmax;
int nw = 501;

double acc_required = 1e-6;

TEST(pade, lorentzian)
{
    double mu=0.2;
    double gamma=0.3;

    // G(iw) for a model with rectangular DOS
    vector< complex<double> > giw(N/2);
    for(int i=0; i<N/2; i++){
        double omega = double(2*i+1) * M_PI / beta;
        giw[i] = 1./ ( IMAG*omega + mu + IMAG*gamma );
    }

    class Gf gf;
    gf.set_Giw(giw, Gf::FERMION, beta);

    // Matsubara frequencies
    vector<double> w(nw), rho(nw);
    for(int i=0; i<nw; i++){
        w[i] = wmin + double(i) * (wmax-wmin) / double(nw-1);
        rho[i] = gf.rho(w[i]);
    }

    // create reference data
    ofstream fout;
    fout.open("pade_lo.dat");
    for(int i=0; i<nw; i++)
        fout << scientific << setprecision(10) << rho[i] << std::endl;
    fout.close();

    // read reference data
    std::vector<double> ref;
    read_file_to_vector("data/ref/pade_lo.dat", ref);

    // compare with reference data
    double diff = diff_vector(rho, ref);
    printf("%.3e\n", diff);

    EXPECT_GT(acc_required, diff);
}


TEST(pade, rectangular)
{
    // G(iw) for a model with rectangular DOS
    vector< complex<double> > giw(N/2);
    for(int i=0; i<N/2; i++){
        double omega = double(2*i+1) * M_PI / beta;
        giw[i] = -IMAG * atan(1./omega);
        // truncate some digits
        // Output of atan slightly depend on environment, which results in failure of make test
        giw[i] = complex<double>( (float)giw[i].real(), (float)giw[i].imag() );
}

    class Gf gf;
    gf.set_Giw(giw, Gf::FERMION, beta);

    // Matsubara frequencies
    vector<double> w(nw), rho(nw);
    for(int i=0; i<nw; i++){
        w[i] = wmin + double(i) * (wmax-wmin) / double(nw-1);
        rho[i] = gf.rho(w[i]);
    }

    // create reference data
    ofstream fout;
    fout.open("pade_rect.dat");
    for(int i=0; i<nw; i++)
        fout << scientific << setprecision(10) << rho[i] << std::endl;
    fout.close();

    // read reference data
    std::vector<double> ref;
    read_file_to_vector("data/ref/pade_rect.dat", ref);

    // compare with reference data
    double diff = diff_vector(rho, ref);
    printf("%.3e\n", diff);

    EXPECT_GT(acc_required, diff);
}
