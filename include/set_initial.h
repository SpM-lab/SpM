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
#ifndef _SET_INITIAL_HEADER
#define _SET_INITIAL_HEADER
#include <string>
#include <map>
#include <vector>
#include "spm_param.h"
#include <algorithm>

class SetInitial {
private:
    struct OmegaInfo {
        int NW;
        double omega_min;
        double omega_max;
    };
    struct FileInfo {
        std::string filein_G;
        int col;
        int output_interval;
        std::string fileout_spec;
        int print_level;
        FileInfo(){
          print_level=1;
          output_interval=1;
        }
    };

    struct CalcInfo {
        std::string statistics;
        double beta;
    };

    int argc;
    std::vector<std::string> vargv;
    SPM_Param param;
    SPM_Flags flags;

    std::map<std::string, std::string> mapForKeyWordToValue;
    std::map<std::string, bool> mapForKeyWordToRead;

    template<class T>
    std::string argv_or_defaultvalue(int n, T value);
    bool SetDefaultValue();
    void SetInputValue();
    void RegisterMap(std::string _keyword, std::string _value);
    void DeleteSpace(std::string &_item){
        size_t c;
        while((c = _item.find_first_of(" \"\'")) != std::string::npos){
            _item.erase(c,1);
        }
    }

    std::string GetValue(std::string _keyword){
        std::transform(_keyword.begin(), _keyword.end(), _keyword.begin(), toupper);
        return mapForKeyWordToValue[_keyword];
    }

public:
    SPM_Param GetParam() { return param; }
    SPM_Flags GetFlags() { return flags; }

    OmegaInfo omegaInfo;
    FileInfo fileInfo;
    CalcInfo calcInfo;
    bool AddDefaulValueMap(char* _filename);
    bool ReadParam(char* _filename);
    bool InputFromArgs(int argc, char *argv[]);
    void PrintInfo();
};


#endif
