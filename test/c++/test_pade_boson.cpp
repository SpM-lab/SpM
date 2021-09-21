#include "gtest.h"
#include "Gf.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iterator>

using namespace std;
const complex<double> IMAG = complex<double>(0, 1.);

bool isodd(int n) { return (n/2)*2 != n;}

template <typename T>
void read_file_to_vector(const char *filename, vector<T> &out)
{
    out.clear();
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
        // std::cerr << i << " " << vec1[i] << " " << vec2[i] << " " << abs(vec1[i]-vec2[i]) << " " << dif_max << std::endl;
    }
    return dif_max/vec1.size();
}

double beta=100.0;

double wmax = 2.0;

double acc_required = 1e-3;

double test_boson()
{
    // read data
    std::vector<double> gtau;
    read_file_to_vector("data/boson/Gtau.in", gtau);
    std::vector<double> ref;
    read_file_to_vector("data/boson/Gtau.in.dos", ref);

    for(int i=0; i<gtau.size(); ++i){
        gtau[i] *= -1.0;
    }

    Gf gf;
    gf.set_Gtau(gtau, Gf::BOSON, beta, 0.0);

    int nw = ref.size();
    vector<double> w(nw), rho(nw);
    for(int i=0; i<nw; i++){
        w[i] = i*wmax /(nw-1);
        rho[i] = gf.rho(w[i]);
        ref[i] *= w[i];
    }

    // compare with reference data
    double diff = diff_vector(rho, ref);
    printf("%.3e\n", diff);

    // EXPECT_GT(acc_required, diff);
    return diff;
}

TEST(pade, pade_boson)
{
    double diff = test_boson();
    EXPECT_GT(acc_required, diff);
}

