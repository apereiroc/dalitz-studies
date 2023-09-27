//
// Created by Asier Pereiro Castro on 27/9/23.
//

#include <iostream>
#include "FitParameters.h"

int main() {

  const double originalVal1 = 0.4;
  const double originalVal2 = -4.9;
  const double originalVal3 = 10.9928;

  const double err = 0.001;

  FitParameters mn_param;

  mn_param.Add("Param1", originalVal1, err);
  mn_param.Add("Param2", originalVal2, err);
  mn_param.Add("Param3", originalVal3, err);

  mn_param.AddBlindingString("Param3", "testBlinding");


  std::cout << "Initial parameters: " << mn_param << "\n\n";

  mn_param.Blind();

  std::cout << "Blinded parameters: " << mn_param << "\n\n";

  mn_param.Unblind();

  std::cout << "Unblinded parameters: " << mn_param << std::endl;

  return 0;
}