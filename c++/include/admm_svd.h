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
/*
 compute sparse solution x of the equation y = A.x by solving
   min_x { (1/2/lambda) ||y-A.x||_2^2 + ||V^t.x||_1 }
     subject to x_i >= 0 & sum_i x_i = 'sum_x'
 where V^t is defined by SVD of A:  A = U.W.V^t

  x: N-dimensional vector
  y: M-dimensional vector
  A: (M * N) matrix
  B: (K * N) matrix, where K is number of singular values retained
*/

#ifndef _ADMM_H
#define _ADMM_H

#include <vector>
#include "cpplapack.h"

class SVD_matrix;
struct admm_result{
	// w/ sv  : results in SV basis
	// w/o sv : results in omega-tau basis
	std::vector<double> x, xsv, z1, z1sv, z2, z2sv, z3, z3sv;
	std::vector<double> y, ysv, y_recovered_x, y_recovered_z1, ysv_recovered_x, ysv_recovered_z1;
};
struct admm_info{
	double res1_pri, res1_dual, res2_pri, res2_dual;  // residual errors
	double res3_pri, res3_dual; // spmpade
	double mse, mse_full;
	int l0_norm;
	double l1_norm, sum_x_calc, negative_weight;
	int iter;
};

class admm_svd{
public:
	admm_svd(CPPL::dgematrix &A, double SV_min=0);  // A
	~admm_svd();

	// [optional]
	void set_sumrule(double sum_x);
	void set_rho_ref(std::vector<double>const &rho_ref, std::vector<double> const & rho_ref_weight);
	void set_omega_coeff(std::vector<double> const& omega_coeff);
	void set_nonnegative(bool _flag);
	void set_fileout_iter(const std::string filename);  // output convergence in a file.  unset if filename=""
	void set_print_level(int);  // 0: none,  1: results,  2: verbose

	// [required]
	void set_coef(double lambda, double penalty1=1.0, double penalty2=1.0, double penalty3=1.0, bool flag_penalty_auto=false);
	void set_y(const std::vector<double> &y);
	int get_svd_num(){return svd_num;}

	// [optional]
	void clear_x();

	// [required]
	int solve(double tolerance=1e-6, int max_iter=1000, int min_iter=10);  // return 0 if converged, and 1 if not converged

	struct admm_result result;
	struct admm_info info;

	double validate(const std::vector<double> &y);  // return MSE

private:
	const CPPL::dgematrix A;
	SVD_matrix *svd;

	bool flag_nonnegative;
	bool flag_sumrule;
	bool flag_rho_ref;
	double sum_x;
	std::vector<double> rho_ref;
	std::vector<double> rho_ref_weight;

	CPPL::dcovector x, Vx, WVx, z1, u1, z2, u2;  // 1 for L1 norm, 2 for non-negativity
	CPPL::dcovector z3, u3; // SpM-Pade
	CPPL::dcovector y, y_sv;

	double regulariz, penalty1, penalty2;
	double penalty3; // SpMPade
	bool flag_penalty_auto;
	static const int PENALTY_UPDATE_INTERVAL;

	std::vector<double> omega_coeff; // omega

	int print_level;
	bool flag_fileout_iter;
	std::string file_iter;

	int svd_num;
	void pre_update();  // This function must be called before calling update_x, and must be recalled when one of values of lambda, penalty1, penalty2 are changed
	void update_x();

	CPPL::dcovector transform_x2sv(CPPL::dcovector const& v);
	CPPL::dcovector transform_sv2x(CPPL::dcovector const& v);

	// quantities used in functions update_x and set_y (set in function pre_update)
	struct quantities_for_update{
		CPPL::dgematrix Y, B;  // diagonal matrix
		CPPL::dgbmatrix Pade_B3;
		CPPL::dgematrix Pade_B0, C;
		CPPL::dcovector Yy, w;
		CPPL::dcovector Pade_y;
		CPPL::drovector v_row;
        CPPL::dcovector rho_w;
		double sum_Vw;
	} pre;
};

#endif
