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
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "cpplapack.h"
#include <cassert>
#include <string>
#include <string.h>
#include "spm_core.h"
#include "kernel.h"
#include "set_initial.h"

using namespace std;

//=============================================================================
// Read one line in the file 'fp' and store the data in 'param[]'.
// An empty line and the line starting with '#' are skipped.
// return: number of data stored in params[]

typedef vector<double> type_gtau;

static int read_data(FILE *fp, double *params) {
  char str[65536];

  do {
    if (fgets(str, 65536, fp) == NULL) return 0;
// 		printf("%s", str);
// 		printf("%d\n", strlen(str));
  } while (str[0] == '#' || strlen(str) == 2);

  int n = 0;
  char *p_str = strtok(str, " ");
  while (p_str != NULL) {
    sscanf(p_str, "%lf", &params[n++]);
    p_str = strtok(NULL, " ");
  }

  return n;
}

static type_gtau read_Gtau(const char *filename, int col) {
  type_gtau Gtau;  // vector<double>

  FILE *fp;
  if ((fp = fopen(filename, "r")) != NULL) {
    double data[128];
    while (read_data(fp, data) > col) {
      Gtau.push_back(data[col]);
    }
    fclose(fp);
  } else {
    printf("file '%s' not found\n", filename);
  }

  return Gtau;
}

static bool if_exist(const char *filename) {
  FILE *fp;
  if ((fp = fopen(filename, "r")) == NULL) return false;
  else return true;
}

//[TODO] Read from File lists
static void read_Gtau_postfix(const char *filename, int col, vector<type_gtau> &g_array) {
  g_array.clear();
  char str[128];
  for (int i = 0; i < 10; i++) {
    sprintf(str, "%s.%01d", filename, i);  // .0, .1, ...
    if (if_exist(str)) g_array.push_back(read_Gtau(str, col));
  }
  for (int i = 0; i < 100; i++) {
    sprintf(str, "%s.%02d", filename, i);  // .01 etc
    if (if_exist(str)) g_array.push_back(read_Gtau(str, col));
  }
  for (int i = 0; i < 1000; i++) {
    sprintf(str, "%s.%03d", filename, i);  // .001 etc
    if (if_exist(str)) g_array.push_back(read_Gtau(str, col));
  }
}

static type_gtau &operator+=(type_gtau &g1, const type_gtau &g2) {
  assert(g1.size() == g2.size());
  for (unsigned i = 0; i < g1.size(); i++) g1[i] += g2[i];
  return g1;
}

static type_gtau &operator/=(type_gtau &g1, double a) {
  for (unsigned i = 0; i < g1.size(); i++) g1[i] /= a;
  return g1;
}

static type_gtau average(vector<type_gtau> &g_array) {
  type_gtau ave(g_array[0].size(), 0.);
  for (unsigned i = 0; i < g_array.size(); i++) ave += g_array[i];
  ave /= double(g_array.size());
  return ave;
}

static type_gtau average_rest(vector<type_gtau> &g_array, unsigned k) {
  type_gtau ave(g_array[0].size(), 0.);
  for (unsigned i = 0; i < g_array.size(); i++) if (i != k) ave += g_array[i];
  ave /= double(g_array.size() - 1);
  return ave;
}

// return -vector
static vector<double> operator-(vector<double> &v) {
  vector<double> u(v.size());
  for (unsigned i = 0; i < v.size(); i++) u[i] = -v[i];
  return u;
}

// matrix-vector product, vec1 = mat*vec2
static vector<double> operator*(CPPL::dgematrix &mat, vector<double> &vec) {
  int m = mat.m;
  int n = mat.n;
  assert(vec.size() == n);
  vector<double> u(m, 0);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) u[i] += mat(i, j) * vec[j];
  }
  return u;
}

//=============================================================================

int main(int argc, char *argv[]) {
  SetInitial setInitial;
  int iret = 0;
  if(argc !=2){
      std::cout << "input file \""<<argv[1]<<"\" does not exist." << std::endl;
      exit(-1);
  }
  else{
      if (setInitial.ReadParam(argv[1]) == false) {
          std::cout << "fail to read input file \""<<argv[1]<<"\"." << std::endl;
          exit(-1);
      }
  }
  setInitial.PrintInfo();

  std::string outputDir = "./output/";
  std::string command = "mkdir -p " + outputDir;
  system(command.c_str());

  // Read Gtau
  printf("\nReading Gtau...\n");
  SPM_Flags flags = setInitial.GetFlags();

  type_gtau Gtau;  // vector<double>
  vector<type_gtau> Gtau_samples;  // for cross validation
  vector<type_gtau> Gtau_rest;  // for cross validation
  if (!flags.validation) {
      Gtau = read_Gtau(setInitial.fileInfo.filein_G.c_str(), setInitial.fileInfo.col);
  } else {
      // cross validation
      // For Obuchi-Kabashima method
      Gtau = read_Gtau(setInitial.fileInfo.filein_G.c_str(), setInitial.fileInfo.col);

      //read_Gtau_postfix(setInitial.fileInfo.filein_G.c_str(), setInitial.fileInfo.col, Gtau_samples);
      //printf("size=%lu\n", Gtau_samples.size());
      //Gtau = average(Gtau_samples);
      //for (unsigned i = 0; i < Gtau_samples.size(); i++) Gtau_rest.push_back(average_rest(Gtau_samples, i));
  }

  type_gtau Gtau_sigma;
  if(setInitial.GetParam().pade.eta != 0.0){
    if(!std::isfinite(setInitial.calcInfo.sigma)){
      if(setInitial.fileInfo.colsigma < 0){
        printf("Neither G_sigma nor column_sigma is set.\n");
        exit(-1);
      }
      Gtau_sigma = read_Gtau(setInitial.fileInfo.filein_Gsigma.c_str(),
                                       setInitial.fileInfo.colsigma);
    }else{
      Gtau_sigma.resize(Gtau.size(), setInitial.calcInfo.sigma);
    }
  }else{
    Gtau_sigma.resize(Gtau.size(), 0.0);
  }

  Kernel kernel;
  int M = Gtau.size();
  std::string statistics = setInitial.calcInfo.statistics;
  double beta = setInitial.calcInfo.beta;
  vector<double> tau = kernel.mesh_linear(0, beta, M);
  printf(" M = %d\n", M);

  if (setInitial.fileInfo.print_level >= 2) {
    FILE *fp = fopen((outputDir + "mesh_tau.dat").c_str(), "w");
    for (unsigned i = 0; i < tau.size(); i++) fprintf(fp, "%d %.5e %.5e\n", i, tau[i], tau[i] / beta);
    fclose(fp);
    printf(" 'mesh_tau.dat'\n");
  }

  //---------------------------------------------------------
  // Set omega mesh
  printf("\nSetting omega...\n");
  vector<double> omega = kernel.mesh_linear(setInitial.omegaInfo.omega_min, setInitial.omegaInfo.omega_max,
                                            setInitial.omegaInfo.NW);
  double d_omega = omega[1] - omega[0];
  printf(" d_omega = %lf\n", d_omega);
  if (setInitial.fileInfo.print_level >= 2) {
    FILE *fp = fopen((outputDir + "mesh_omega.dat").c_str(), "w");
    for (unsigned i = 0; i < omega.size(); i++) fprintf(fp, "%d %.5e %.5e\n", i, omega[i], omega[i] * beta);
    fclose(fp);
    printf(" 'mesh_omega.dat'\n");
  }

  //---------------------------------------------------------

  SPM_Core spm_core;
  SPM_Param spm_param = setInitial.GetParam();
  //Get lambda
  vector<double> lambda;
  lambda = kernel.mesh_log(spm_param.lambda.lbegin, spm_param.lambda.lend, spm_param.lambda.Nl);

  //Set parameters
  spm_core.SetParameters(setInitial.GetParam());
  spm_core.SetFlags(setInitial.GetFlags());

  printf("\nConstructing kernel matrix...\n");
  vector<vector<double> > A;
  vector<double> vmse, vmse_full, l1_norm, valid;
  vector<int> l0_norm;
  int errinfo = 0;

  //TODO: Make Kernel should be treated as a private function
  errinfo = kernel.MakeKernelLinear(statistics, beta, tau, omega, A, M, setInitial.omegaInfo.NW);
  if (errinfo != 0) {
    return errinfo;
  }

  printf("\nSolving the equation...\n");
  int l_valid;
  errinfo = spm_core.SolveEquation(statistics, beta, A, Gtau, Gtau_sigma, lambda, omega);
  if (errinfo != 0) {
    return errinfo;
  }

  spm_core.GetLambdaOpt(lambda, &l_valid);

  vector<double> spec;
  spm_core.GetSpectrum(spec);
  if (statistics == "boson"){
      for(int i = 0; i < omega.size(); i++){
          double x = exp(fabs(omega[i])*beta);
          double factor = (1.0 - x)/ (1.0 + x);
          if(fabs(omega[i]) < 1e-15) factor = beta/2.0;
          else{
              if (omega[i] > 0) factor *= -1.0;
              factor /= omega[i];
          }
          spec[i] *= factor;
      }
  }
  spm_core.GetResults(vmse, vmse_full, l0_norm, l1_norm, valid);
  //---------------------------------------------------------
  // Save spectrum
  {
    FILE *fp = fopen((outputDir + setInitial.fileInfo.fileout_spec).c_str(), "w");
    fprintf(fp, "# lambda=%.3e  (l=%d)\n", lambda[l_valid], l_valid);
    if (statistics == "fermion") {
        for (int i = 0; i < setInitial.omegaInfo.NW; i++) {
            fprintf(fp, "%.5e %.5e\n", omega[i], spec[i] / d_omega);
        }
    }
    else{
        for (int i = 0; i < setInitial.omegaInfo.NW; i++) {
            fprintf(fp, "%.5e %.5e %.5e\n", omega[i], spec[i]*omega[i] / d_omega, spec[i] / d_omega);
        }
    }
    fclose(fp);
    printf("'%s'\n", setInitial.fileInfo.fileout_spec.c_str());
  }

  //---------------------------------------------------------
  // Save lambda dependence
  {
    string fileout_lambda("lambda_dep.dat");
    FILE *fp = fopen((outputDir + fileout_lambda).c_str(), "w");
    for (unsigned l = 0; l < lambda.size(); l++) {
      fprintf(fp, "%.5e %.5e %.5e %d %.5e %.5e\n", lambda[l], vmse[l], vmse_full[l], l0_norm[l], l1_norm[l], valid[l]);
    }
    fclose(fp);
    printf("'%s'\n", fileout_lambda.c_str());
  }
}
