//
// Created by Asier Pereiro Castro on 13/9/23.
//

#pragma once

#include <string>
#include <sstream>

#define MIN_BLIND_VALUE 1000.0
#define MAX_BLIND_VALUE 20000.0

// Function to convert a hash value to a double
//double hashValueToDouble(size_t hashValue) {
//  // Use a stringstream to convert the hash value to a string
//  std::stringstream ss;
//  ss << hashValue;
//
//  // Convert the string to a double
//  double result;
//  ss >> result;
//  return result;
//}
//
//// Function to hash a std::string and convert it to a double
//double hashStringToDouble(const std::string &inputString) {
//  std::hash<std::string> hasher;
//  size_t hashValue = hasher(inputString);
//  return hashValueToDouble(hashValue);
//}

double hashValueToDouble(size_t hashValue, double minValue, double maxValue) {
  // Calculate the range size and map the hash value to the desired range
  double range = maxValue - minValue;
  double scaledValue = minValue + (static_cast<double>(hashValue) /
                                   static_cast<double>(std::numeric_limits<size_t>::max())) *
                                  range;
  return scaledValue;
}

// Function to hash a std::string, reduce the range, and convert it to a double
double hashStringToDouble(const std::string &inputString, double minValue,
                          double maxValue) {
  std::hash<std::string> hasher;
  size_t hashValue = hasher(inputString);
  return hashValueToDouble(hashValue, minValue, maxValue);
}


class Blinder {
public:
  static double getOffset(const std::string &blindingString) {
    return hashStringToDouble(blindingString, MIN_BLIND_VALUE, MAX_BLIND_VALUE);
  }
};
