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
#include <stdio.h>
#include <fstream> // for debug
#include <spm_core.h>
#include "common.h"
#include <Gf.h>
#include <cassert>
#include <cmath>

#include <random>

double SPM_Core::integrate(std::vector<double> &y, double width) {
  long unsigned int n = y.size();
  double s = (y[0] + y[n - 1]) / 2.;
  for (long unsigned int i = 1; i < n - 1; i++) s += y[i];
  return s * width / double(n - 1);
}

int SPM_Core::SetKernel(std::vector<std::vector<double> > &_AIn) {

  unsigned int ANrow = _AIn.size();
  unsigned int ANcolumn = _AIn[0].size();

  if (ANrow < 1 || ANcolumn < 1) {
    return ERROR_INVALID_MATRIX_SIZE;
  }
  A.resize(ANrow, ANcolumn);
  for (unsigned int i = 0; i < ANrow; i++) {
    for (unsigned int j = 0; j < ANcolumn; j++) {
      A(i, j) = _AIn[i][j];
    }
  }
  return ERROR_SUCCESS;
}

int SPM_Core::SetParametersLambda(int _NLambda, double _lbegin, double _lend) {
  //TODO: Check value
  param.lambda.Nl = _NLambda;
  param.lambda.lbegin = _lbegin;
  param.lambda.lend = _lend;
  return ERROR_SUCCESS;
}

int SPM_Core::SetParametersAdmm(double _penalty, double _tolerance, int _max_iter) {
  //TODO: Check value
  param.admm.penalty = _penalty;
  param.admm.tolerance = _tolerance;
  param.admm.max_iter = _max_iter;
  return ERROR_SUCCESS;
}

int SPM_Core::SetFlagValidation(bool _flag_validation) {
  //TODO: Check value
  flags.validation = _flag_validation;
  return ERROR_SUCCESS;
}

int SPM_Core::SetParametersSVD(double _SV_min) {
  //TODO: Check value
  param.svd.sv_min = _SV_min;
  return ERROR_SUCCESS;
}

int SPM_Core::SolveEquation(
        std::string _StatisticsType,
        double _Beta,
        std::vector<std::vector<double> > &_AIn,
        std::vector<double> &_Gtau,
        std::vector<double> &_Gtau_error,
        std::vector<double> &_lambda,
        std::vector<double> &_omega) {

  this->StatisticsType = _StatisticsType;
  this->beta = _Beta;
  std::vector<double> y = _Gtau;  // input
  if (this->StatisticsType == "fermion") {
    for (unsigned int i = 0; i < y.size(); i++) {
      y[i] *= -1.0;
    }
  }
  std::vector<double> _omega_coeff(_omega);
  double sum_G;
  if (this->StatisticsType == "fermion") {
    sum_G = -_Gtau[0] - _Gtau[_Gtau.size() - 1];
    std::fill(_omega_coeff.begin(), _omega_coeff.end(), 1.0);
  } else if (_StatisticsType == "boson") {
    sum_G = _Gtau[0] + _Gtau[_Gtau.size() - 1];
    const size_t nw = _omega_coeff.size();
    for(size_t iw = 0; iw < nw; iw++){
      const double bw = _Beta * _omega[iw];
      if(bw > 0.0){
        _omega_coeff[iw] = -std::expm1(-bw) / (1.0 + std::exp(-bw));
      }else{
        _omega_coeff[iw] = std::expm1(bw) / (1.0 + std::exp(bw));
      }
    }
  } else {
    std::cerr << "Error: StatisticsType must be fermion or boson." << std::endl;
    return -1;
  }

  std::vector<double> rho_ref(_omega.size());
  std::vector<double> rho_ref_weight(_omega.size());
  std::vector<double> rhopade(_omega.size());

  // Pade from G(tau)
  {
    std::vector<double> gt2(_Gtau.begin(), _Gtau.end()-1);
    Gf gf;
    if(_StatisticsType=="fermion"){
      // The tail (4th argument) is not used in SpM
      gf.set_Gtau(gt2, Gf::FERMION, _Beta, 1.0);
    }else{
      // The tail (4th argument) is not used in SpM
      for(int i=0; i<gt2.size(); ++i){
        gt2[i] *= -1.0;
      }
      gf.set_Gtau(gt2, Gf::BOSON, _Beta, 0.0);
    }

    for(int i=0; i<_omega.size(); ++i){
      rhopade[i] = gf.rho(_omega[i]);
    }
  }

  if(param.pade.eta != 0.0){
    // Pade from G(tau) + Gaussian noise
    std::mt19937 rng(137);
    std::normal_distribution<> randn(0.0,1.0);
    const int Nsample = param.pade.nsample;
    for(int isample=0; isample<Nsample; ++isample){
      std::vector<double> gt2(_Gtau.begin(), _Gtau.end()-1);
      const int Nt = gt2.size();
      for(int i=0; i<Nt; ++i){
        gt2[i] += _Gtau_error[i]*randn(rng);
      }
      Gf gf;
      if(_StatisticsType=="fermion"){
        // The tail (4th argument) is not used in SpM
        gf.set_Gtau(gt2, Gf::FERMION, _Beta, 1.0);
      }else{
        for(int i=0; i<gt2.size(); ++i){
          gt2[i] *= -1.0;
        }
        // The tail (4th argument) is not used in SpM
        gf.set_Gtau(gt2, Gf::BOSON, _Beta, 0.0);
      }
      for(int i=0; i<_omega.size(); ++i){
        double r = gf.rho(_omega[i]);
        rho_ref[i] += r;
        rho_ref_weight[i] += r*r;
      }
    }
  }

  const double dw = _omega[1] - _omega[0];

  FILE *fp = fopen(("./output/" + param.pade.filename).c_str(), "w");
  if(param.pade.eta == 0.0){
    for(int i=0; i<_omega.size(); ++i){
      std::fprintf(fp, "%.15lf %.15lf\n", _omega[i], rhopade[i]);
    }
  }else{
    const int Nsample = param.pade.nsample;
    for(int i=0; i<_omega.size(); ++i){
      rho_ref_weight[i] -= rho_ref[i]*rho_ref[i]/Nsample;
      rho_ref_weight[i] /= (Nsample-1); 
      rho_ref[i] /= Nsample;
      auto sigma = std::sqrt(rho_ref_weight[i]);

      if(param.pade.eta < 0.0){
        // constant weight
        rho_ref_weight[i] = -param.pade.eta;
      }else{
        auto r = std::abs(sigma / rho_ref[i]);
        rho_ref_weight[i] = param.pade.eta/(1.0+r*r);
      }
      std::fprintf(fp, "%.15lf %.15lf %.15lf %.15lf %.15lf\n",
             _omega[i],
             rhopade[i],
             rho_ref[i],
             sigma,
             rho_ref_weight[i]
             );
      rho_ref[i] *= dw;
    }
  }
  fclose(fp);

  return SolveEquationCore(_AIn, y, _omega, _omega_coeff, _lambda, sum_G, rho_ref, rho_ref_weight);
}

int SPM_Core::SolveEquationCore(
        std::vector<std::vector<double> > &_AIn,
        std::vector<double> &_y,
        std::vector<double> &_omega,
        std::vector<double> &_omega_coeff,
        std::vector<double> &_lambda,
        const double _sum_G,
        std::vector<double> &_rho_ref,
        std::vector<double> &_rho_ref_weight
        ) {
  SetKernel(_AIn);
  // Solve the equations
  lambda = _lambda;
  std::vector<double> y = _y;  // input

  bool flag_penalty_auto = param.admm.penalty < 0 ? true : false;
  param.admm.penalty = std::abs(param.admm.penalty);

  admm_svd D(A, param.svd.sv_min);
  D.set_nonnegative(flags.nonnegative);

  //D.set_print_level(print_level);
  std::string outputOrgDir = "./output/";

  for (unsigned int l = 0; l < lambda.size(); l++) {
    printf("\n==================================================\n");
    printf(" lambda = %.3e\n", lambda[l]);

    std::string prefix;  // directory name
    {  // set prefix and creat directory
      char str[64];
      std::string outputDir = outputOrgDir + "lambda/";
      sprintf(str, "lambda_%.2e/", lambda[l]);
      // printf("%s\n", str);
      prefix = outputDir + str;
      std::string command = "mkdir -p " + prefix;
      system(command.c_str());
    }
    D.set_fileout_iter(prefix + "iter.dat");

    //core
    D.set_sumrule(_sum_G);

    D.set_rho_ref(_rho_ref, _rho_ref_weight);
    D.set_omega_coeff(_omega_coeff);
    D.set_coef(lambda[l], param.admm.penalty, param.admm.penalty, param.admm.penalty, flag_penalty_auto);
    D.set_y(y);
    D.solve(param.admm.tolerance, param.admm.max_iter);

    omega = _omega;
    double d_omega = omega[1] - omega[0];
    print_results_admm(D.result, omega, d_omega, prefix);

    result.push_back(D.result);
    info.push_back(D.info);

    // cross validation
    /*
    if( flags.validation){
      int K=Gtau_samples.size();
      double diff=0;
      for(int k=0; k<K; k++){
        printf("\n----- k=%d\n", k);
        D.set_y(-Gtau_samples[k]);
        D.solve(admm_tolerance, admm_max_iter);
        diff += D.validate(-Gtau_rest[k]);
      }
      diff /= (double)K;
      printf("valid=%.5e\n", diff);
      valid.push_back(diff);
    }
    */
  }
  printf("\n-----\n");

  // determine lambda
  param.lambda.lvalid = 0;
  if (flags.validation) {
      std::vector<double> looc(param.lambda.Nl);
      int M = D.get_svd_num();
      char filename[64] = "find_lambda_opt_CV.dat";
      FILE *fp = fopen((outputOrgDir + filename).c_str(), "w");
      fprintf(fp, "# log10(lambda)  log10(E_LOCC) log10(lambda-E_LOCC) log10(E_LOCC-lambda) M S(lambda)\n");
      for (int l = 0; l < param.lambda.Nl; l++){
          int Slambda=info[l].l0_norm;
          double x = M/static_cast<double>(M-Slambda);
          looc[l] = info[l].mse*x*x;
          fprintf(fp, " %.15e %.15e %.15e %.15e %d %d\n",
                  log10(lambda[l]), log10(looc[l]),
                  log10(lambda[l]-looc[l]), log10(looc[l]-lambda[l]),
                  M, Slambda);
      }
      fclose(fp);
      printf("'%s'\n", filename);

      // get position of the minimum element
      param.lambda.lvalid = std::distance(looc.begin(), std::min_element(looc.begin(), looc.end()));
      valid.resize(lambda.size());
      fill(valid.begin(), valid.end(), 0.);
  } else {
    std::vector<double> mse, diff, log_f;
    for (int l = 0; l < param.lambda.Nl; l++) mse.push_back(info[l].mse);
    param.lambda.lvalid = find_kink(lambda, mse, diff, log_f);
    if (param.lambda.lvalid != -1) {
      // print
      char filename[64] = "find_lambda_opt.dat";
      FILE *fp = fopen((outputOrgDir + filename).c_str(), "w");
      fprintf(fp, "# log(x)  diff  log(y)  log(f(x))  [all in log10 scale]\n");
      for (int l = 0; l < param.lambda.Nl; l++) {
        fprintf(fp, " %.5e %.5e %5e %5e\n", log10(lambda[l]), diff[l] / log(10), log10(mse[l]), log_f[l] / log(10));
      }
      fclose(fp);
      printf("'%s'\n", filename);
    } else {
      param.lambda.lvalid = 0;
    }

    valid.resize(lambda.size());
    fill(valid.begin(), valid.end(), 0.);
  }
  printf("*** optimal lambda = %.3e  (l=%d)\n", lambda[param.lambda.lvalid], param.lambda.lvalid);

  // copy optimal lambda folder
  std::string prefix_opt;  // directory name
  {  // set prefix and creat directory
    char str[64];
    std::string orgDir = outputOrgDir + "lambda/";
    std::string toDir = outputOrgDir + "lambda_opt";
    sprintf(str, "lambda_%.2e/", lambda[param.lambda.lvalid]);
    prefix_opt = orgDir + str;
    std::string command = "cp -rf " + prefix_opt + " " + toDir;
    system(command.c_str());
  }

  return 0;
}

int SPM_Core::find_kink(std::vector<double> &x, std::vector<double> &y, std::vector<double> &diff,
                        std::vector<double> &log_f) {
  if (y.size() < 3) return -1;

  // f = a*x^b
  // log(f) = b*log(x) + log(a)
  double y1 = y.front();
  double y2 = y.back();
  double x1 = x.front();
  double x2 = x.back();
  double b = log(y1 / y2) / log(x1 / x2);
  double log_a = log(y2) - b * log(x2);

  // diff = log(f) - log(y)
  diff.assign(y.size(), 0);
  log_f.assign(y.size(), 0);
  for (unsigned l = 0; l < y.size(); l++) {
    log_f[l] = (b * log(x[l]) + log_a);
    diff[l] = (log_f[l] - log(y[l]));
  }
  // find maximum
  int l_opt = distance(diff.begin(), max_element(diff.begin(), diff.end()));  // get position of the maximum element
  return l_opt;
}

void SPM_Core::GetSpectrum(std::vector<double> &_spectrum) {
  if(flags.nonnegative==true || flags.refrho ) {
    _spectrum = result[param.lambda.lvalid].z2;
  }
  else{
    _spectrum = result[param.lambda.lvalid].z1;
  }
}

void SPM_Core::GetLambdaOpt(std::vector<double> &_lambda, int *_opt_l) {
  _lambda = lambda;
  *_opt_l = param.lambda.lvalid;
}

void SPM_Core::GetResults(std::vector<double> &_vmse, std::vector<double> &_vmse_full,
                          std::vector<int> &_vl0_norm, std::vector<double> &_vl1_norm,
                          std::vector<double> &_valid) {
  _vmse.resize(info.size());
  _vmse_full.resize(info.size());
  _vl0_norm.resize(info.size());
  _vl1_norm.resize(info.size());
  _valid.resize(info.size());
  for (unsigned i = 0; i < info.size(); i++) {
    _vmse[i] = info[i].mse;
    _vmse_full[i] = info[i].mse_full;
    _vl0_norm[i] = info[i].l0_norm;
    _vl1_norm[i] = info[i].l1_norm;
    _valid[i] = valid[i];
  }
}

void SPM_Core::print_results_admm(admm_result &r, std::vector<double> w, double dw, std::string prefix) {
  FILE *fp;
  std::string filename;

  // x in SV basis
  filename = prefix + "x_sv.dat";
  fp = fopen(filename.c_str(), "w");
  for (unsigned i = 0; i < r.xsv.size(); i++) {
    fprintf(fp, "%d %.5e %.5e %.5e %.5e\n", i, r.xsv[i], r.z1sv[i], r.z2sv[i], r.z3sv[i]);
    // fprintf(fp, "%d %.5e %.5e %.5e\n", i, r.xsv[i]/dw, r.z1sv[i]/dw, r.z2sv[i]/dw);
    // fprintf(fp, "%d %.5e %.5e %.5e\n", i, r.xsv[i]/sqrt(dw), r.z1sv[i]/sqrt(dw), r.z2sv[i]/sqrt(dw));
  }
  fclose(fp);
  printf("'%s'\n", filename.c_str());

  // x in tau-omega basis
  filename = prefix + "x_tw.dat";
  fp = fopen(filename.c_str(), "w");
  for (unsigned i = 0; i < r.x.size(); i++) {
    fprintf(fp, "%.5e %.5e %.5e %.5e %.5e\n", w[i], r.x[i] / dw, r.z1[i] / dw, r.z2[i] / dw, r.z3[i] / dw);
  }
  fclose(fp);
  printf("'%s'\n", filename.c_str());

  // y in SV basis
  filename = prefix + "y_sv.dat";
  fp = fopen(filename.c_str(), "w");
  for (unsigned i = 0; i < r.ysv.size(); i++) {
    fprintf(fp, "%d %.5e %.5e %.5e\n", i, r.ysv[i], r.ysv_recovered_x[i], r.ysv_recovered_z1[i]);
  }
  fclose(fp);
  printf("'%s'\n", filename.c_str());
  // y in tau-omega basis
  filename = prefix + "y_tw.dat";
  fp = fopen(filename.c_str(), "w");
  for (unsigned i = 0; i < r.y.size(); i++) {
    fprintf(fp, "%.5e %.5e %.5e %.5e\n", static_cast<double>(i) / (r.y.size()-1), r.y[i], r.y_recovered_x[i],
            r.y_recovered_z1[i]);
  }
  fclose(fp);
  printf("'%s'\n", filename.c_str());
}
