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
#include <cassert>
#include "admm_svd.h"
// included <vector>, "cpplapack.h"
#include "common.h"
#include "SVD_matrix.h"

using namespace std;
// const string admm_svd::file_iter = "admm_iter.dat";
const int admm_svd::PENALTY_UPDATE_INTERVAL = 20;

//=============================================================================
// class admm_svd

///
/// \param A_in [in]
/// \param SV_min [in]
admm_svd::admm_svd(CPPL::dgematrix &A_in, double SV_min) : A(A_in), flag_sumrule(false), flag_nonnegative(false),
                                                           flag_fileout_iter(false), print_level(1) {
  svd = new SVD_matrix(A, SV_min);
  if (print_level >= 0) {
    std::string outputDir = "./output/";
    std::string command = "mkdir -p " + outputDir;
    system(command.c_str());
    svd->OutputSVD(outputDir + "SV.dat");
    if (print_level >= 2)
      svd->print_basis(outputDir + "SVD_U.dat", outputDir + "SVD_V.dat", 10);
  }

  int N = A.n;  // size of x
  int K = svd->get_rank();  // # of SV bases

  x.resize(K);
  z1.resize(K);
  u1.resize(K);
  z2.resize(N);
  u2.resize(N);
  clear_x();
}

admm_svd::~admm_svd() {
  delete svd;
}

void admm_svd::set_sumrule(double sum_x) {
  this->sum_x = sum_x;
  flag_sumrule = true;
}

void admm_svd::set_nonnegative() {
  flag_nonnegative = true;
}

void admm_svd::set_fileout_iter(string filename) {
  if (filename.empty()) {
    flag_fileout_iter = false;
  } else {
    flag_fileout_iter = true;
    file_iter = filename;
  }
}

void admm_svd::set_print_level(int print_level) {
  this->print_level = print_level;
}
//=============================================================================

void admm_svd::set_coef(double regulariz, double penalty1, double penalty2, bool flag_penalty_auto) {
  this->regulariz = regulariz;
  this->penalty1 = penalty1;
  this->penalty2 = flag_nonnegative ? penalty2 : 0;
  this->flag_penalty_auto = flag_penalty_auto;

  pre_update();
}

void admm_svd::pre_update() {
  int K = svd->get_rank();

  CPPL::dgbmatrix diag(K, K, 0, 0);  // diagonal matrix  D^{-1}
  for (int i = 0; i < K; i++) diag(i, i) = 1. / (penalty1 + penalty2 + pow(svd->S(i, i), 2) / regulariz);

  pre.Y = diag * svd->S / regulariz;
  pre.B = diag * penalty1;
  pre.C = diag * svd->VT * penalty2;

  CPPL::dcovector p = dcovector_all1(A.n);
  pre.w = diag * svd->VT * p;  // D^{-1} * V^t * p
  pre.v_row = CPPL::t(p) * CPPL::t(svd->VT);  // p^t * V
  pre.sum_Vw = pre.v_row * pre.w;  // < V*w > = p^t * V * w

  {  // check matrix size
    int M = A.m, N = A.n;
    assert(pre.Y.m == K && pre.Y.n == K);  // Y: K*K
    assert(pre.B.m == K && pre.B.n == K);  // B: K*K
    assert(pre.C.m == K && pre.C.n == N);  // C: K*N
    assert(pre.w.l == K);
    assert(pre.v_row.l == K);
  }
}

void admm_svd::set_y(const vector<double> &y_in) {
  assert(A.m == y_in.size());  // y: M
  y = vec2cppl_col(y_in);
  y_sv = svd->transform_y2sv(y);  // y' = U^t * y
  pre.Yy = pre.Y * y_sv;  // Y * y'
}

void admm_svd::clear_x() {
  x.zero();
  z1.zero();
  u1.zero();
  z2.zero();
  u2.zero();
}

//=============================================================================

// soft thresholding function S_a(y)
static double soft_threshold(double a, double y) {
  assert(a > 0);
  if (fabs(y) <= a) return 0;  // -a <= y <= a
  if (y > 0) return y - a;  // y > a
  else return y + a;  // y < -a
}

static CPPL::dcovector soft_threshold(double a, CPPL::dcovector y) {
  CPPL::dcovector r(y.l);
  for (int i = 0; i < y.l; i++) {
    r(i) = soft_threshold(a, y(i));
  }
  return r;
}

void admm_svd::update_x() {
  // update x
  x = pre.Yy;
  x += pre.B * (z1 - u1);
  if (flag_nonnegative) x += pre.C * (z2 - u2);

  if (flag_sumrule) {
    double nu = (sum_x - pre.v_row * x) / pre.sum_Vw;  // pre.v_row*x　=　<V*x>
    x += nu * pre.w;
  }
  Vx = svd->transform_sv2x(x);  // V x'

  // update z1 and u1
  z1 = soft_threshold(1. / penalty1, x + u1);
  u1 += x - z1;

  // update z2 and u2
  if (flag_nonnegative) {
    z2 = positive(Vx + u2);  // projection
    u2 += Vx - z2;
  }
}

// update penalty and vector u
// r: primary residual,  s: dual residual
// Ref: Eq.(3.13) in Boyd et al., Foundations and Trends in Machine Learning Vol. 3, No. 1 (2010) 1–122
static bool update_penalty(double r, double s, double &penalty, CPPL::dcovector &u) {
  const double MAX_RESIDUAL_RATIO = 10;
  const double PENALTY_INCREMENT = 2;
  double fac = PENALTY_INCREMENT;
  if (r > s * MAX_RESIDUAL_RATIO) fac = fac;  // increase penalty if primary residual is too large
  else if (s > r * MAX_RESIDUAL_RATIO) fac = 1. / fac;  // decrease penalty if dual residual is too large
  else return false;  // no change

  penalty *= fac;
  u /= fac;
  return true;
}

// input: primary residual error per element
static bool if_converge(double r, double tol) {
  if (r < tol) return true;
  return false;
}

int admm_svd::solve(double tolerance, int max_iter, int min_iter) {
  // iteration
  int iter = 0, counter = 0;
  double diff_x = 0;
  info.res1_pri = info.res1_dual = info.res2_pri = info.res2_dual = 0;
  bool flag_converge = false;
  int update_interval1 = PENALTY_UPDATE_INTERVAL;
  int update_interval2 = PENALTY_UPDATE_INTERVAL;
  FILE *fp = NULL;
  if (flag_fileout_iter) {
    fp = fopen(file_iter.c_str(), "w");
    fprintf(fp, "# iter  diff(x,x_old)  res1_pri res1_dual res2_pri res2_dual  RMSE L1_norm sum(x) negative_weight\n");
  }
  do {
    CPPL::dcovector x_old(x);  // copy
    CPPL::dcovector z1_old(z1);
    CPPL::dcovector z2_old(z2);

    //---------------------------------------------------------------------
    // UPDATE x, z1, u1, z2, u2, Vx
    update_x();

    //---------------------------------------------------------------------
    // ERRORS AND COST FUNCTIONS

    // residual errors
    info.res1_pri = norm_l2(x - z1);  // propto sqrt(K)
    info.res1_dual = norm_l2(z1 - z1_old) * penalty1;  // propto sqrt(N)
    if (flag_nonnegative) {
      info.res2_pri = norm_l2(Vx - z2);  // propto sqrt(N)
      info.res2_dual = norm_l2(z2 - z2_old) * penalty2;  // propto sqrt(N)  // ****check
    }

    // quantities to check
    diff_x = norm_l1(x - x_old);
    info.mse = norm_l2_sq(y_sv - svd->S * x);  // MSE computed in SV basis
    info.l1_norm = norm_l1(x);
    // sum_x_calc = sum_vector(Vx);
    info.sum_x_calc = pre.v_row * x;  // < V*x >
    info.negative_weight = sum_vector(positive(-Vx));  // < P_+(V.x) >

    //---------------------------------------------------------------------
    // PRINT OUT
    if (print_level >= 2) {
      if (flag_nonnegative && if_negative(Vx)) printf("*");  // non-negativity is not satisfied
      else printf(".");  // OK
      fflush(stdout);
    }

    if (flag_fileout_iter && iter) {
      fprintf(fp, "%d %.6e", iter, diff_x);
      fprintf(fp, " %.6e %.6e %.6e %.6e", info.res1_pri, info.res1_dual, info.res2_pri, info.res2_dual);
      fprintf(fp, " %.6e %.6e %.6e %.6e", info.mse, info.l1_norm, info.sum_x_calc, info.negative_weight);
      fprintf(fp, "\n");
    }

    //---------------------------------------------------------------------
    // UPDATE PENALTY
    if (flag_penalty_auto) {
      counter++;
      // printf("examine penalty\n");
      if (counter % update_interval1 == 0) {
        bool flag_update = false;
        flag_update |= update_penalty(info.res1_pri, info.res1_dual, penalty1, u1);
        if (flag_nonnegative) {
          flag_update |= update_penalty(info.res2_pri, info.res2_dual, penalty2, u2);
        }
        if (flag_update) {
          // set_matrix_SVD();
          pre_update();
          if (print_level >= 2) {
            printf("\n update penalty1=%lf  penalty2=%lf  (iter=%d)\n", penalty1, penalty2, iter + 1);
          }
          // fprintf(fp, "\n");
          counter = 0;
          update_interval1 *= 2;  // extend interval
        }
      }
    }

    //---------------------------------------------------------------------
    // CHECK CONVERGENCE
    // criteria: absolute tolerance (primary residual error per element)
    // Ref: Eq.(3.12) in Boyd et al., Foundations and Trends in Machine Learning Vol. 3, No. 1 (2010) 1–122
    flag_converge = if_converge(info.res1_pri / sqrt(z1.l), tolerance);
    flag_converge &= if_converge(info.res2_pri / sqrt(x.l), tolerance);

    //---------------------------------------------------------------------
    iter++;
  } while ((!flag_converge && iter < max_iter) || iter < min_iter);
  fclose(fp);

  info.mse_full = norm_l2_sq(y - A * Vx);
  info.iter = iter;

  if (print_level) {
    printf(" iteration finish\n");
    if (flag_converge) printf(" converged\n");
    else printf(" not converged\n");

    if (flag_fileout_iter) printf(" | '%s'\n", file_iter.c_str());
    printf(" | iter            = %d\n", info.iter);
    printf(" | diff(x, x_old)  = %.3e\n", diff_x);
    printf(" |--- residual errors\n");
    printf(" | L1 regul. (pri) = %.3e\n", info.res1_pri);
    printf(" |          (dual) = %.3e\n", info.res1_dual);
    printf(" | x>=0      (pri) = %.3e\n", info.res2_pri);
    printf(" |          (dual) = %.3e\n", info.res2_dual);
    printf(" |--- target quantities\n");
    printf(" | MSE(y-A.x)      = %.3e (computed in SV basis)\n", info.mse);
    printf(" | MSE(y-A.x)      = %.3e (computed in tau-omega basis)\n", info.mse_full);
    printf(" | || x' ||_1      = %.3e\n", info.l1_norm);
    printf(" | sum(x)          = %.6e", info.sum_x_calc);
    if (flag_sumrule) printf("  (error = %.2e)", fabs(sum_x - info.sum_x_calc));
    printf("\n");
    printf(" | negative weight = %.3e\n", info.negative_weight);
    printf(" |--- penalty parameters\n");
    printf(" |   for L1 regul. = %.3lf\n", penalty1);
    printf(" |   for x>=0      = %.3lf\n", penalty2);
    printf(" |---\n");
  }

  {  // copy results
    result.x = cppl2vec(Vx);
    result.xsv = cppl2vec(x);
    result.z1 = cppl2vec(svd->transform_sv2x(z1));  // V * z1
    result.z1sv = cppl2vec(z1);
    result.z2 = cppl2vec(z2);
    result.z2sv = cppl2vec(svd->transform_x2sv(z2));  // V^t * z2

    result.y = cppl2vec(y);
    result.ysv = cppl2vec(y_sv);
    result.y_recovered_x = cppl2vec(svd->transform_sv2y(svd->S * x));  // U y'
    result.y_recovered_z1 = cppl2vec(svd->transform_sv2y(svd->S * z1));  // U y'
    result.ysv_recovered_x = cppl2vec(svd->S * x);  // y'
    result.ysv_recovered_z1 = cppl2vec(svd->S * z1);  // y'
  }

  return flag_converge ? 0 : 1;
}

double admm_svd::validate(const vector<double> &y_test) {
  CPPL::dcovector y2 = vec2cppl_col(y_test);
  CPPL::dcovector y2_sv = svd->transform_y2sv(y2);
  // return norm_l2_sq(y2 - A*Vx);
  return norm_l2_sq(y2_sv - svd->S * x);
  // return norm_l2_sq(y2 - A*svd->transform_sv2x(z1));
  // return norm_l2_sq(y2_sv - svd->S*z1);
}
//=============================================================================
