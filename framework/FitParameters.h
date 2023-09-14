#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"
#include "Blinder.h"

using json = nlohmann::ordered_json;


/**
   Enhance MnUserParameters functionality
 */
class FitParameters : public ROOT::Minuit2::MnUserParameters {
private:
  std::vector<std::string> latexNames;
  std::vector<std::pair<unsigned int, double>> blindingOffsets; // First: index, Second: blinding offset

  void UpdatePar(const ROOT::Minuit2::FunctionMinimum &min) {
    const ROOT::Minuit2::MnUserParameters &par_min = min.UserParameters();

    for (unsigned int i = 0; i < par_min.Params().size(); ++i) {
      this->SetValue(i, par_min.Value(i));
      this->SetError(i, par_min.Error(i));
      if (par_min.Parameter(i).IsFixed())
        this->Fix(i);
    }
  }

public:
  inline void LoadParFromJSON(const std::string &file) {
    std::ifstream f(file);

    json parameters = json::parse(f);

    latexNames.reserve(parameters.size());

    for (auto &[p, data]: parameters.items()) {
      const std::string name = p;
      const double val = data["val"];
      const bool fixed = data["fixed"];
      const std::string latex = data["latex"];

      if (!data["err"].is_null())
        this->Add(name, val, data["err"]);
      else
        this->Add(name, val, 99999.0);
      if (fixed)
        this->Fix(name);

      if (!data["ll"].is_null() and !data["ul"].is_null()) {
        this->SetLimits(name, data["ll"], data["ul"]);
      } else if (!data["ll"].is_null()) {
        this->SetLowerLimit(name, data["ll"]);
      } else if (!data["ul"].is_null()) {
        this->SetUpperLimit(name, data["ul"]);
      }

      latexNames.push_back(latex);
    }
  }

  void SaveParToJSON(const std::string &file) const {

    json parameters;

    for (std::size_t i = 0; i < this->Params().size(); i++) {

      const std::string name = this->Parameter(i).GetName();
      json p;

      p["val"] = this->Parameter(i).Value();
      p["err"] = this->Parameter(i).Error();
      p["fixed"] = this->Parameter(i).IsFixed();

      if (this->Parameter(i).HasLowerLimit())
        p["ll"] = this->Parameter(i).LowerLimit();
      else
        p["ll"] = nullptr;

      if (this->Parameter(i).HasUpperLimit())
        p["ul"] = this->Parameter(i).UpperLimit();
      else
        p["ul"] = nullptr;

      p["latex"] = latexNames.at(i);

      parameters[name] = p;
    }

    std::ofstream output(file);

    output << std::setw(4) << parameters << std::endl;
  }

  void LoadBlindingScheme(const std::string &file) {
    std::ifstream f(file);

    json json_file = json::parse(f);

    const unsigned int npars = this->Params().size();

    if (npars == 0)
      std::cerr
              << "No parameters found. Please add parameters before adding their blinding scheme\n";

    for (auto &[paramName, data]: json_file.items()) {
      const std::string &blindingString = data["blinding_string"];

      const unsigned int idx = this->Index(paramName);

      if (idx > npars)
        std::cerr << "Parameter not found\n";

      const double offset = Blinder::getOffset(blindingString);

      blindingOffsets.push_back({idx, offset});
    }
  }

  void SaveParToJSON(const std::string &file,
                     const ROOT::Minuit2::FunctionMinimum &min) {
    if (!min.IsValid()) {
      std::cout << "WARNING: Note that minimum is invalid" << std::endl;
    }
    UpdatePar(min);
    SaveParToJSON(file);
  }

  void Blind() {
    for (auto &blindingOffset: blindingOffsets) {
      const unsigned int idx = blindingOffset.first;
      const double offset = blindingOffset.second;
      const double currentVal = this->Value(idx);

      this->SetValue(idx, currentVal + offset);
    }
  }

  void Unblind() {
    for (auto &blindingOffset: blindingOffsets) {
      const unsigned int idx = blindingOffset.first;
      const double offset = blindingOffset.second;
      const double currentVal = this->Value(idx);

      this->SetValue(idx, currentVal - offset);
    }
  }

  void FixAllParams() {
    for (std::size_t i = 0; i < this->Params().size(); i++) {
      this->Fix(i);
    }
  }

  FitParameters &operator+=(const FitParameters &other) {
    for (unsigned int i = 0; i < other.Params().size(); ++i) {
      const std::string &name = other.GetName(i);
      const double value = other.Value(i);
      const double error = other.Error(i);

      this->Add(name, value, error);

      double ll = 0.0, ul = 0.0;
      if (other.Parameter(i).HasLowerLimit())
        ll = other.Parameter(i).LowerLimit();
      if (other.Parameter(i).HasUpperLimit())
        ul = other.Parameter(i).UpperLimit();

      if (other.Parameter(i).HasLowerLimit() &&
          other.Parameter(i).HasUpperLimit())
        this->SetLimits(name, ll, ul);
      else if (other.Parameter(i).HasLowerLimit() &&
               !other.Parameter(i).HasUpperLimit())
        this->SetLowerLimit(name, ll);
      else if (!other.Parameter(i).HasLowerLimit() &&
               other.Parameter(i).HasUpperLimit())
        this->SetUpperLimit(name, ul);
      else {
      }

      if (other.Parameter(i).IsFixed())
        this->Fix(other.Name(i));
    }

    return *this;
  }

  void SetAndFix(const std::string &name, const double &value) {
    this->SetValue(name, value);
    this->Fix(name);
  }

  void SetAndFix(const unsigned int &number, const double &value) {
    this->SetValue(number, value);
    this->Fix(number);
  }
};

FitParameters operator+(FitParameters lhs, const FitParameters &rhs) {
  return lhs += rhs;
}
