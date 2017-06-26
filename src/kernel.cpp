/* SPM  -  Sparse Modeling tool */
/* Copyright (C) 2017 Junya Otsuki, Kazuyoshi Yoshimi, Hiroshi Shinaoka, Masayuki Ohzeki*/

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
#include "kernel.h"
#include <cassert>
#include "errorcode.h"

using namespace std;

///
/// \param _Beta
/// \param _Tau
/// \param _Omega
/// \param _Aout
/// \param _ANrow
/// \param _ANcolumn
/// \return
int Kernel::MakeKernelLinear
        (
                const std::string _StatisticsType,
                const double _Beta,
                const vector<double> &_Tau,
                const vector<double> &_Omega,
                vector<vector<double> > &_Aout,
                const int _ANrow,
                const int _ANcolumn
        ) {
  if (_Tau.size() == 0 || _Omega.size() == 0) {
    return (ERROR_INVALID_SIZE);
  }
  omega = _Omega;
  if (_ANrow < 1 || _ANcolumn < 1) {
    return (ERROR_INVALID_MATRIX_SIZE);
  }
  _Aout.resize(_ANrow);
  for (int i = 0; i < _ANrow; i++) {
    _Aout[i].resize(_ANcolumn);
  }

  double (*kernel)(double, double, double) = _StatisticsType == "fermion" ? Kernel::kernel_f : Kernel::kernel_b;

  for (int j = 0; j < _ANcolumn; j++) {
    for (int i = 0; i < _ANrow; i++) {
      _Aout[i][j] = kernel(_Beta, _Tau[i], omega[j]);
    }
  }
  return ERROR_SUCCESS;
}

///
/// \param x_min
/// \param x_max
/// \param N
/// \return
vector<double> Kernel::mesh_linear(double x_min, double x_max, int N) {
  if (N == 1) return vector<double>(1, x_min);
  vector<double> x(N);
  for (int i = 0; i < N; i++) x[i] = x_min + double(i) * (x_max - x_min) / double(N - 1);
  return x;
}

///
/// \param x_min
/// \param x_max
/// \param N
/// \return
vector<double> Kernel::mesh_log(double x_min, double x_max, int N) {
  if (N == 1) return vector<double>(1, x_min);
  assert(x_min > 0);
  assert(x_max > 0);
  vector<double> x_log = mesh_linear(log(x_min), log(x_max), N);
  vector<double> x(N);
  for (int i = 0; i < N; i++) x[i] = exp(x_log[i]);
  return x;
}

double Kernel::kernel_f(double beta, double tau, double omega) {
  double bw = beta * omega;
  double tw = tau * omega;
  if (omega >= 0) return exp(-tw) / (1. + exp(-bw));
  else return exp(bw - tw) / (1. + exp(bw));
}

double Kernel::kernel_b(double beta, double tau, double omega) {
  double bw = beta * omega;
  double tw = tau * omega;
  if (fabs(bw) < 1e-6) return exp(-tw) / beta / (1. - bw / 2. + bw * bw / 6. - bw * bw * bw / 24.);
  else if (omega > 0) return exp(-tw) / (1. - exp(-bw)) * omega;
  else return exp(bw - tw) / (exp(bw) - 1.) * omega;
}

