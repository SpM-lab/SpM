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

#ifndef _SVD_MATRIX_H
#define _SVD_MATRIX_H

#include <vector>
#include "cpplapack.h"

class SVD_matrix {
    friend class admm;
    friend class admm_svd;

public:
    SVD_matrix(CPPL::dgematrix A, double SV_min = 0);  // copy the matrix A because it is destroyed in SVD
    int get_rank();

    CPPL::dgematrix At();  // A^t
    CPPL::dgematrix At_A();  // A^t * A
    CPPL::dgematrix At_A_inv();  // (A^t * A)^{-1}
    CPPL::dgematrix At_A_inv_At();  // (A^t * A)^{-1} * A^t
    void rearrange_col(std::vector<int> &);  // rearrange column of A (actaully VT)
    void print_basis(std::string _file_U, std::string _file_V, int);

    CPPL::dcovector transform_x2sv(const CPPL::dcovector &);  // V^t * x  (to SV basis)
    CPPL::dcovector transform_y2sv(const CPPL::dcovector &);  // U^t * y  (to SV basis)
    CPPL::dcovector transform_sv2x(const CPPL::dcovector &);  // V * x'  (to original basis)
    CPPL::dcovector transform_sv2y(const CPPL::dcovector &);  // U * y'  (to original basis)
    void OutputSVD(std::string _file_SVD);
private:
    CPPL::dcovector S_temp;
    CPPL::dgematrix U_temp, VT_temp;
    // K <= M,N
    CPPL::dgbmatrix S, S_inv;  // singular values, K*K diagonal matrix
    CPPL::dgematrix U, VT, VT_full;  // M*K, K*N
    std::string file_SVD;
    std::string file_U, file_V;
};

#endif
