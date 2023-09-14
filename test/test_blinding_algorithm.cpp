//
// Created by Asier Pereiro Castro on 14/9/23.
//

#include <iostream>
#include "Blinder.h"

int main() {

  const std::vector<std::string> inputs = {
          "",
          "a",
          "b",
          "ab",
          "Hello",
          "Hello world!",
          "aaaaaaaaaaaaaaaaa",
          "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
          "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
  };

  for (const auto &input: inputs) {
    const double output = Blinder::getOffset(input);
    std::cout << "Output for {" << input << "}: " << output << std::endl;
  }

  return 0;
}