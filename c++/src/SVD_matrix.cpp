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
#include "SVD_matrix.h"
#include "common.h"
#include <cassert>

using namespace std;

/*
const string SVD_matrix::file_SVD = "SV.dat";
const string SVD_matrix::file_U = "SVD_U.dat";
const string SVD_matrix::file_V = "SVD_V.dat";
*/


void SVD_matrix::OutputSVD(std::string _file_SVD) {
	this->file_SVD = _file_SVD;
	std::cout << this->file_SVD << std::endl;
	FILE *fp = fopen(file_SVD.c_str(), "w");
	for (int i = 0; i < S_temp.l; i++) fprintf(fp, "%d %.5e\n", i, S_temp(i));
	fclose(fp);
	printf(" | '%s'\n", file_SVD.c_str());

}

///
/// \param A
/// \param SV_min
SVD_matrix::SVD_matrix(CPPL::dgematrix A, double SV_min) {
	printf("# SVD\n");

	// SVD
	A.dgesvd(S_temp, U_temp, VT_temp);  // A is destroyed

	// low-rank approximation
	int K = min(A.m, A.n);
	printf(" | largest SV = %.2e\n", S_temp(0));
	for (int k = 0; k < min(A.m, A.n); k++) {
		// printf("%d %.3e\n", k, S(k));
		if (S_temp(k) < SV_min) {
			K = k;
			break;
		}
	}
	printf(" | # of SVs retained = %d   ( SV > %.2e )\n", K, SV_min);
	S.resize(K, K, 0, 0);  // diagonal matrix
	S_inv.resize(K, K, 0, 0);
	for (int i = 0; i < K; i++) {
		S(i, i) = S_temp(i);
		S_inv(i, i) = 1. / S_temp(i);
	}
	U = low_rank_matrix(U_temp, A.m, K);
	VT = low_rank_matrix(VT_temp, K, A.n);
	VT_full = VT;
}

int SVD_matrix::get_rank() {
	return S.m;
}

CPPL::dgematrix SVD_matrix::At() {
	return CPPL::t(VT) * S * CPPL::t(U);
}

CPPL::dgematrix SVD_matrix::At_A() {
	return CPPL::t(VT) * (S * S) * VT;
}

CPPL::dgematrix SVD_matrix::At_A_inv() {
	return CPPL::t(VT) * (S_inv * S_inv) * VT;
}

CPPL::dgematrix SVD_matrix::At_A_inv_At() {
	return CPPL::t(VT) * S_inv * CPPL::t(U);
}

void SVD_matrix::rearrange_col(vector<int> &index) {
	VT = arrange_matrix_col(VT_full, index);
}

// print first N bases of U and V matrices
void SVD_matrix::print_basis(std::string _file_U, std::string _file_V, int N) {
	this->file_U = _file_U;
	this->file_V = _file_V;
	{  // U matrix
		FILE *fp = fopen(file_U.c_str(), "w");
		int m = U.m;
		if (U.n < N) {
			N = U.n;
		}
		for (int i = 0; i < m; i++) {
			fprintf(fp, "%d", i);
			for (int j = 0; j < N; j++) fprintf(fp, " %.5e", U(i, j));
			// CPPL::dcovector u(U.col(j));
			// for(int i=0; i<m; i++)  fprintf(fp, "%d %.5e\n", i, u(i));
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" | '%s'\n", file_U.c_str());
	}
	{  // V matrix
		FILE *fp = fopen(file_V.c_str(), "w");
		int m = VT.n;
		for (int i = 0; i < m; i++) {
			fprintf(fp, "%d", i);
			for (int j = 0; j < N; j++) fprintf(fp, " %.5e", VT(j, i));
			// CPPL::dcovector u(U.col(j));
			// for(int i=0; i<m; i++)  fprintf(fp, "%d %.5e\n", i, u(i));
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf(" | '%s'\n", file_V.c_str());
	}
}

CPPL::dcovector SVD_matrix::transform_x2sv(const CPPL::dcovector &x) {
	assert(VT.n == x.l);
	return VT * x;  // V^t * x
}

CPPL::dcovector SVD_matrix::transform_y2sv(const CPPL::dcovector &y) {
	assert(U.m == y.l);
	return CPPL::t(U) * y;  // U^t * y
}

CPPL::dcovector SVD_matrix::transform_sv2x(const CPPL::dcovector &xsv) {
	assert(VT.m == xsv.l);
	return CPPL::t(VT) * xsv;  // V * x'
}

CPPL::dcovector SVD_matrix::transform_sv2y(const CPPL::dcovector &ysv) {
	assert(U.n == ysv.l);
	return U * ysv;  // U * y'
}

