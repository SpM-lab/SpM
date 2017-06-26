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
#ifndef _SPM_PARAM_HEADER
#define _SPM_PARAM_HEADER

class SPM_Param{
 private:
	struct Lambda{
		int Nl;
		double lbegin;
		double lend;
		int lvalid;
    double dlambda;
		Lambda(){
			Nl=1;
			lbegin = 1e-1;
			lend = 1e+0;
			lvalid=0;
			dlambda=-1;
		}
	};

	struct Admm{
		double penalty;
		double tolerance;
		int max_iter;
		bool flag_penalty_auto;
		Admm(){
			penalty=10.;
			tolerance = 1e-6;
			max_iter=1000;
			flag_penalty_auto=false;
		}
	};

	struct SVD{
		double sv_min;
		SVD(){
			sv_min=0;
		}
	};
	
 public:
    Lambda lambda;
    Admm admm;
    SVD svd;
    
};

struct SPM_Flags{
  bool validation;
};
  
#endif
