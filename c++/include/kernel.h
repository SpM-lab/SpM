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

#ifndef _KERNEL_H
#define _KERNEL_H

#include <string>
#include <vector>
#include "cpplapack.h"

class Kernel{

private:
	CPPL::dgematrix A;
	std::vector <double> omega;
	std::vector <double> lambda;
	std::vector <double> valid;
		static double kernel_f(double beta, double tau, double omega);
		static double kernel_b(double beta, double tau, double omega);

public:

	std::vector<double> mesh_linear(double x_min, double x_max, int N);
	std::vector<double> mesh_log(double x_min, double x_max, int N);

		int MakeKernelLinear(
					const std::string _StatisticsType,
					const double _Beta,
					const std::vector<double> &_Tau,
					const std::vector <double> &_Omega,
					std::vector<std::vector<double> > &_Aout,
					const int _ANrow, const int _ANcolum
	);

};

#endif
