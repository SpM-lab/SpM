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
 perform SVD, A = U*S*V^t, and store those matrices
 SVs smaller than 'SV_min' are dropped.
 All SVs are retained if SV_min=0.
*/

#ifndef _SPM_CORE_H
#define _SPM_CORE_H

#include <vector>
#include "cpplapack.h"
#include "errorcode.h"
#include "admm_svd.h"
#include "spm_param.h"

class SPM_Core {

private:
		CPPL::dgematrix A;
		SPM_Param param;
		SPM_Flags flags;
		std::vector<double> lambda;
		std::vector<double> valid;
		std::vector<double> omega;
		std::vector<admm_result> result;
		std::vector<admm_info> info;

		double integrate(std::vector<double> &y, double width);

		int find_kink(std::vector<double> &x, std::vector<double> &y, std::vector<double> &diff, std::vector<double> &log_f);

		void print_results_admm(admm_result &r, std::vector<double> w, double dw, std::string prefix);

public:

		int SetKernel(std::vector<std::vector<double> > &_AIn);

		void SetParameters(const SPM_Param &_param){param=_param;}

		void SetFlags (const SPM_Flags &_flags){flags=_flags;}

		int SetParametersSVD(double _SVD_min);

		int SetParametersLambda(int _NLambda, double _lbegin, double _lend);

		int SetParametersAdmm(double _penalty, double _tolerance, int _max_iter);

		int SetFlagValidation(bool _flag_validation);

		void GetLambdaOpt(std::vector<double> &_lambda, int *_opt_l);

		void GetSpectrum(std::vector<double> &_spectrum);

		void GetResults(std::vector<double> &_vmse, std::vector<double> &_vmse_full, std::vector<double> &_vl1_norm,
										std::vector<double> &_valid);

		int SolveEquation(
						std::string _StatisticsType,
						double _Beta,
						std::vector<std::vector<double> > &_AIn,
						std::vector<double> &_Gtau,
            std::vector<double> &_lambda,
            std::vector<double> &_omega);

		int SolveEquationCore(
						std::vector<std::vector<double> > &_AIn,
						std::vector<double> &_Gtau,
						std::vector<double> &_omega,
            std::vector<double> &_lambda,
						const double _sum_G
		);
};

#endif
