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
#ifndef _COMMON_H
#define _COMMON_H

#include <vector>
#include "cpplapack.h"

// converter
CPPL::dcovector vec2cppl_col(const std::vector<double> &);  // vector -> CPPL
std::vector<double> cppl2vec(const CPPL::dcovector &);  // CPPL -> vector
std::vector<double> cppl2vec(const CPPL::drovector &);  // CPPL -> vector

int norm_l0(const CPPL::dcovector &);  // L0 norm
double norm_l1(const CPPL::dcovector &);  // L1 norm
double norm_l2(const CPPL::dcovector &);  // L2 norm
double norm_l2_sq(const CPPL::dcovector &);  // square of L2 norm
double sum_vector(const CPPL::dcovector &);  // sum of vector elements
double innerproduct(const CPPL::dcovector &, const CPPL::dcovector &);  //innerproduct

bool if_negative(const CPPL::dcovector &);  // true if at least one component is x_i<0
bool if_negative(const CPPL::dgematrix &);

CPPL::dcovector dcovector_all1(int n);  // return a covector with elements all being 1
CPPL::drovector drovector_all1(int n);  // return a rovector with elements all being 1

CPPL::dcovector positive(const CPPL::dcovector &);  // set negative elements at zero
CPPL::dgematrix positive(const CPPL::dgematrix &);

// return shrinked matrix
CPPL::dgematrix low_rank_matrix(CPPL::dgematrix &mat, int m, int n);

// return matrix in which the columns are rearranged according to 'index_col' (it can be smaller than the original column-length)
CPPL::dgematrix arrange_matrix_col(const CPPL::dgematrix &mat_full, const std::vector<int> &index_col);

#endif
