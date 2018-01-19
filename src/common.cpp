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
#include <cassert>
#include "common.h"

using namespace std;

// converter: vector -> CPPL
CPPL::dcovector vec2cppl_col(const vector<double> &v) {
	CPPL::dcovector u(v.size());
	for (unsigned i = 0; i < v.size(); i++) {
		u(i) = v[i];
	}
	return u;
}

// converter: CPPL -> vector
vector<double> cppl2vec(const CPPL::dcovector &u) {
	vector<double> v(u.l);
	for (int i = 0; i < u.l; i++) {
		v[i] = u(i);
	}
	return v;
}

vector<double> cppl2vec(const CPPL::drovector &u) {
	vector<double> v(u.l);
	for (int i = 0; i < u.l; i++) {
		v[i] = u(i);
	}
	return v;
}

// L0 norm
int norm_l0(const CPPL::dcovector &v) {
	int s = 0;
	for (int i = 0; i < v.l; i++) {
		if(fabs(v(i)) > pow(10.0, -15)) s++;
	}
	return s;
}



// L1 norm
double norm_l1(const CPPL::dcovector &v) {
	double s = 0;
	for (int i = 0; i < v.l; i++) s += fabs(v(i));
	return s;
}

// L2 norm
double norm_l2(const CPPL::dcovector &v) {
	double s = 0;
	for (int i = 0; i < v.l; i++) s += v(i) * v(i);
	return sqrt(s);
}

// L2 norm ^2
double norm_l2_sq(const CPPL::dcovector &v) {
	double s = 0;
	for (int i = 0; i < v.l; i++) s += v(i) * v(i);
	return s;
}

// sum of vector elements
double sum_vector(const CPPL::dcovector &v) {
	double s = 0;
	for (int i = 0; i < v.l; i++) s += v(i);
	return s;
}

// innerproduct
double innerproduct(const CPPL::dcovector &v1, const CPPL::dcovector &v2) {
	assert(v1.l == v2.l);
	double s = 0;
	for (int i = 0; i < v1.l; i++) s += v1(i) * v2(i);
	return s;
}

// true if at least one component is x_i<0
bool if_negative(const CPPL::dcovector &v) {
	for (int i = 0; i < v.l; i++) {
		if (v(i) < 0) return true;
		// if( v(i)<0 && fabs(v(i))>1e-10 )  return true;  // neglect rounding errors
	}
	return false;
}

bool if_negative(const CPPL::dgematrix &mat) {
	for (int i = 0; i < mat.m; i++) {
		for (int j = 0; j < mat.n; j++) {
			if (mat(i, j) < 0) return true;
			// if( mat(i,j)<0 && fabs(mat(i,j))>1e-10 )  return true;  // neglect rounding errors
		}
	}
	return false;
}

CPPL::dcovector dcovector_all1(int n) {
	CPPL::dcovector v(n);
	for (int i = 0; i < v.l; i++) v(i) = 1.;
	return v;
}

CPPL::drovector drovector_all1(int n) {
	CPPL::drovector v(n);
	for (int i = 0; i < v.l; i++) v(i) = 1.;
	return v;
}

CPPL::dcovector positive(const CPPL::dcovector &v) {
	CPPL::dcovector u(v.l);
	for (int i = 0; i < v.l; i++) u(i) = v(i) > 0 ? v(i) : 0;
	return u;
}

CPPL::dgematrix positive(const CPPL::dgematrix &A) {
	CPPL::dgematrix B(A.m, A.n);
	for (int i = 0; i < A.m; i++) {
		for (int j = 0; j < A.n; j++) B(i, j) = A(i, j) > 0 ? A(i, j) : 0;
	}
	return B;
}

// return shrinked matrix
CPPL::dgematrix low_rank_matrix(CPPL::dgematrix &mat, int m, int n) {
	assert(m <= mat.m);
	assert(n <= mat.n);
	CPPL::dgematrix r(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) r(i, j) = mat(i, j);
	}
	return r;
}

// return matrix in which the columns are rearranged according to 'index_col' (it can be smaller than the original column-length)
CPPL::dgematrix arrange_matrix_col(const CPPL::dgematrix &mat_full, const vector<int> &index_col) {
	assert(mat_full.n >= index_col.size());
	CPPL::dgematrix mat(mat_full.m, index_col.size());
	for (int i = 0; i < mat.m; i++) {
		for (int j = 0; j < mat.n; j++) {
			mat(i, j) = mat_full(i, index_col[j]);
		}
	}
	return mat;
}
