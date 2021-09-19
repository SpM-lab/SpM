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

#include "set_initial.h"
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>


void SetInitial::PrintInfo() {
  printf("\nParameters:\n");
  printf(" Statistics = %s\n", calcInfo.statistics.c_str());
  printf(" beta = %lf\n", calcInfo.beta);
  //printf(" filein_G = '%s' \n", fileInfo.filein_G.c_str(), fileInfo.col);
  printf(" filein_G = '%s' \n", fileInfo.filein_G.c_str());
  printf(" filein_Gsigma = '%s' \n", fileInfo.filein_Gsigma.c_str());
  printf(" fileout_spec = '%s'\n", fileInfo.fileout_spec.c_str());
  printf(" fileout_pade = '%s'\n", param.pade.filename.c_str());

  printf(" omega = [ %lf : %lf ] (%d points)\n", omegaInfo.omega_min, omegaInfo.omega_max, omegaInfo.NW);

  printf("\n ADMM:\n");
  printf(" lambda = [ %.2e : %.2e ] (%d points)\n", param.lambda.lbegin, param.lambda.lend, param.lambda.Nl);
  if (param.lambda.dlambda > 0) {
    printf(" dlambda = %lf \n", param.lambda.dlambda);
  }
  printf(" penalty = %lf\n", param.admm.penalty);
  printf(" penalty_auto = %s\n", param.admm.flag_penalty_auto ? "true" : "false");
  // printf(" regularization = %lf\n", admm_regular);
  printf(" tolerance = %.2e\n", param.admm.tolerance);
  printf(" max_iter = %d\n", param.admm.max_iter);
  printf(" print_level = %d\n", fileInfo.print_level);
  printf(" OutputInterval = %d\n", fileInfo.output_interval);
  printf(" FlagNonNegative = %s\n", flags.nonnegative ? "ON": "OFF");

  if (flags.refrho){
    printf("\n SpM-Pade:\n");
    printf(" PadeEta = %lf\n", param.pade.eta);
    printf(" NSamplePade = %d\n", param.pade.nsample);
  }

  if (flags.validation) {
    printf("\n CROSS VALIDATION\n");

    /*
    // string command = "ls | grep '" + filein_G + ".[0-9]'";
    // cout << command << endl;
    // system(command.c_str());
    //
    // for(int i=0; i<1000; i++){
    // 	char str[128];
    // 	sprintf(str, "%s.%03d", filein_G.c_str(), i);
    // 	printf("%s\n", str);
     // }
     */
  }
}

bool SetInitial::AddDefaulValueMap(char *_filename) {
  std::ifstream fin;
  fin.open(_filename, std::ios::in);
  if (fin.fail()) {
    std::cerr << "fail to open file." << std::endl;
    return false;
  }
  std::string reading_line_buffer;
  std::string item;
  while (!fin.eof()) {
    std::vector<std::string> vreadline;
    std::getline(fin, reading_line_buffer);
    std::stringstream ss(reading_line_buffer);
    while (std::getline(ss, item, ' ') && !item.empty()) {
      vreadline.push_back(item);
    }
    if (vreadline.size() != 2) {
      std::cerr << "Wrong file format." << std::endl;
      return false;
    }
    RegisterMap(vreadline[0], vreadline[1]);
  }
  return true;
}

void SetInitial::RegisterMap(std::string _keyword, std::string _value) {
  DeleteSpace(_keyword);
  DeleteSpace(_value);
  std::transform(_keyword.begin(), _keyword.end(), _keyword.begin(), toupper);
  mapForKeyWordToValue.insert(std::map<std::string, std::string>::value_type(_keyword, _value));
  mapForKeyWordToRead.insert(std::map<std::string,bool>::value_type(_keyword, false));
}


bool SetInitial::SetDefaultValue() {
  RegisterMap("Statistics", "");
  RegisterMap("beta", "");
  RegisterMap("G_sigma", "Inf");
  RegisterMap("column", "");
  RegisterMap("column_sigma", "-1");
  RegisterMap("filein_G", "");
  RegisterMap("filein_G_sigma", "");
  RegisterMap("fileout_spec", "spectrum.out");
  RegisterMap("fileout_pade", "pade.out");
  RegisterMap("NOmega", "");
  RegisterMap("OmegaMin", "");
  RegisterMap("OmegaMax", "");
  RegisterMap("NLambda", "0");
  RegisterMap("LambdaLogBegin", "");
  RegisterMap("LambdaLogEnd", "");
  RegisterMap("LambdaLogMesh", "0");
  RegisterMap("LambdaValid", "0");
  RegisterMap("Penalty", "1.0");
  RegisterMap("Tolerance", "1e-6");
  RegisterMap("MaxIteration", "1000");
  RegisterMap("SVMin", "1e-10");
  RegisterMap("PrintLevel", "1");
  RegisterMap("CrossValidation", "n");
  RegisterMap("OutputInterval", "1");
  RegisterMap("FlagNonNegative", "ON");
  RegisterMap("PadeEta", "0.0");
  RegisterMap("NSamplePade", "30");

  /*
  std::map<std::string, std::string>::iterator it;
  for(it=mapForKeyWordToValue.begin(); it !=mapForKeyWordToValue.end(); it++){
    std::cout<<"Debug: key="<<it->first<< ", value="<< it->second<<std::endl;
  }
  */

  return true;
}

bool SetInitial::ReadParam(char *_filename) {
  SetDefaultValue();
  std::ifstream fin;
  fin.open(_filename, std::ios::in);
  if (fin.fail()) {
    std::cerr << "fail to open file." << std::endl;
    return false;
  }
  std::string reading_line_buffer;
  std::string item;

  while (!fin.eof()) {
    std::vector<std::string> vreadline;
    std::getline(fin, reading_line_buffer);
    std::stringstream ss(reading_line_buffer);
    if (ss.str() != "") {
      while (std::getline(ss, item, '=') && !item.empty()) {
        DeleteSpace(item);
        vreadline.push_back(item);
      }
      if (vreadline[0].size() != 0 && vreadline[0].substr(0, 1) != "#") {
        std::transform(vreadline[0].begin(), vreadline[0].end(), vreadline[0].begin(), toupper);
        if (mapForKeyWordToValue.find(vreadline[0]) == mapForKeyWordToValue.end()) {
          std::cerr << "Wrong keyword: " << vreadline[0] << " in your input file." << std::endl;
          std::map<std::string, std::string>::iterator it;
          std::cerr << "Keyword List" << std::endl;
          std::cerr << "*****************************" << std::endl;
          for (it = mapForKeyWordToValue.begin(); it != mapForKeyWordToValue.end(); it++) {
            std::cerr << "" << it->first << std::endl;
          }
          std::cerr << "*****************************" << std::endl;
          return false;
        }
        //std::cout<<"Debug: "<<vreadline[0]<<",  "<<vreadline[1]<<std::endl;
        mapForKeyWordToValue[vreadline[0]] = vreadline[1];
        mapForKeyWordToRead[vreadline[0]]=true;
      }
    }
  }
  SetInputValue();
  return true;
}

void SetInitial::SetInputValue() {

  { //Check the esseintial keywords
    std::vector<std::string> EssentialKeyWord = {"Statistics", "beta", "column", "filein_G",
                                                 "OmegaMin", "OmegaMax", "LambdaLogBegin", "LambdaLogEnd"};
    for (std::vector<std::string>::iterator it = EssentialKeyWord.begin(); it != EssentialKeyWord.end(); it++) {
      std::transform(it->begin(), it->end(), it->begin(), toupper);
    }

    int iret = 0;
    for (unsigned i = 0; i < EssentialKeyWord.size(); i++) {
      if (mapForKeyWordToRead[EssentialKeyWord[i]] == false) {
        std::cerr << "Error:  " << EssentialKeyWord[i]
                  << " is an essential parameter. Please set the keyword and value in the input file." << std::endl;
        iret = -1;
      }
    }
    if (iret == -1) {
      exit(-1);
    }
  }

  calcInfo.statistics = GetValue("Statistics");
  transform(calcInfo.statistics.begin(), calcInfo.statistics.end(), calcInfo.statistics.begin(), tolower);
  calcInfo.beta = std::stod(GetValue("beta"));
  calcInfo.sigma = std::stod(GetValue("G_sigma"));
  fileInfo.col = std::stoi(GetValue("column"));
  fileInfo.colsigma = std::stoi(GetValue("column_sigma"));
  fileInfo.filein_G = GetValue("filein_G");
  fileInfo.fileout_spec = GetValue("fileout_spec");
  fileInfo.output_interval = std::stoi(GetValue("OutputInterval"));

  {
    std::string buf = GetValue("filein_G_sigma");
    if(buf==""){fileInfo.filein_Gsigma = fileInfo.filein_G;}
    else{fileInfo.filein_Gsigma = buf;}
  }

  omegaInfo.NW = std::stoi(GetValue("NOmega"));
  omegaInfo.omega_min = std::stod(GetValue("OmegaMin"));
  omegaInfo.omega_max = std::stod(GetValue("OmegaMax"));

  param.lambda.Nl = std::stoi(GetValue("NLambda"));
  param.lambda.lbegin = pow(10.0, std::stod(GetValue("LambdaLogBegin")));
  param.lambda.lend = pow(10.0, std::stod(GetValue("LambdaLogEnd")));
  param.lambda.lvalid = std::stoi(GetValue("LambdaValid"));
  param.lambda.dlambda = std::stod(GetValue("LambdaLogMesh"));
  param.admm.penalty = std::stod(GetValue("Penalty"));
  param.admm.tolerance = std::stod(GetValue("Tolerance"));
  param.admm.max_iter = std::stoi(GetValue("MaxIteration"));
  param.svd.sv_min = std::stod(GetValue("SVMin"));

  fileInfo.print_level = std::stoi(GetValue("PrintLevel"));

  flags.validation = (GetValue("CrossValidation") == "y");

  std::string flag=GetValue("FlagNonNegative");
  std::transform(flag.begin(), flag.end(),flag.begin(), toupper);
  if( flag== "ON" or flag=="OFF"){
    flags.nonnegative= (flag == "ON");
  }
  else{
    std::cerr<<"Error:  FlagNonNegative must be ON or OFF."<<std::endl;
    exit(-1);
  }
  /* TODO: Check
  double tmp_lambda=param.lambda.lbegin;
  if(param.lambda.lbegin < param.lambda.lend){
    param.lambda.lbegin=param.lambda.lend;
    param.lambda.lend=tmp_lambda;
  }
  */

  param.pade.filename = GetValue("fileout_pade");
  param.pade.nsample = std::stoi(GetValue("NSamplePade"));
  param.pade.eta = std::stod(GetValue("PadeEta"));
  flags.refrho = param.pade.eta != 0.0;

  if (param.lambda.Nl == 0) {
    if (fabs(param.lambda.dlambda) < pow(10.0, -12)) {
      param.lambda.dlambda = 0.2;
    }
    if (param.lambda.dlambda > 0) {
      param.lambda.Nl = static_cast<int>((log10(param.lambda.lbegin / param.lambda.lend)) / param.lambda.dlambda + 0.5) + 1;
      if (param.lambda.Nl < 1) {
        printf("\ndLambda must be lesser than log10(LambdaBegin/LambdaEnd)=%lf\n",
               (log10(param.lambda.lbegin / param.lambda.lend)));
        exit(-1);
      }
    } else {
      printf("\ndLambda must be greater than 0.");
      exit(-1);
    }
  } else {//NL is defined
    if (fabs(param.lambda.dlambda) > pow(10.0, -12)) {
      printf("\ndLambda and NLambda should not be defined at the same time.");
      exit(-1);
    }
  }
}
