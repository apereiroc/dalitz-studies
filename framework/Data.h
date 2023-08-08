#pragma once

#include "Event.h"
#include "Selection.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>

void LoadToy(const std::string &filename, const std::string &treename,
             std::vector<Event> &data, const bool &has_cache = false) {
  // Read ROOT file
  TFile *const file = new TFile(filename.c_str());
  if (file->IsZombie()) {
    std::cout << "Error opening TFile " << file->GetName() << std::endl;
    exit(1);
  }

  TTree *const tree = dynamic_cast<TTree *>(file->Get(treename.c_str()));
  if (tree == nullptr) {
    std::cout << "There's no tree called " << treename << " in " << filename
              << std::endl;
    file->Close();
    std::exit(1);
  }

  const std::size_t n = tree->GetEntries();

  // Avoid data reallocations
  data.reserve(n);

  // Declare branches
  double time;
  int qtag;

  TLorentzVector *p4Kp = NULL, *p4pim = NULL, *p4Km = NULL, *p4pip = NULL;

  double B_VV_S, B_VV_P, B_VV_D;
  double B_VV_S_CP, B_VV_P_CP, B_VV_D_CP;
  std::array<double, nVVConf> B_VV, B_VV_CP;

  std::complex<double> *S_VV_S = NULL, *S_VV_P = NULL, *S_VV_D = NULL;
  std::complex<double> *S_VV_S_CP = NULL, *S_VV_P_CP = NULL, *S_VV_D_CP = NULL;
  std::array<std::complex<double>, nVVConf> S_VV, S_VV_CP;

  double pdf_gen = 0.0;

  // Set branch addresses
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("qtag", &qtag);

  tree->SetBranchAddress("p4Kp", &p4Kp);
  tree->SetBranchAddress("p4pim", &p4pim);
  tree->SetBranchAddress("p4Km", &p4Km);
  tree->SetBranchAddress("p4pip", &p4pip);

  tree->SetBranchAddress("B_VV_S", &B_VV_S);
  tree->SetBranchAddress("B_VV_P", &B_VV_P);
  tree->SetBranchAddress("B_VV_D", &B_VV_D);
  tree->SetBranchAddress("B_VV_S_CP", &B_VV_S_CP);
  tree->SetBranchAddress("B_VV_P_CP", &B_VV_P_CP);
  tree->SetBranchAddress("B_VV_D_CP", &B_VV_D_CP);

  tree->SetBranchAddress("S_VV_S", &S_VV_S);
  tree->SetBranchAddress("S_VV_P", &S_VV_P);
  tree->SetBranchAddress("S_VV_D", &S_VV_D);
  tree->SetBranchAddress("S_VV_S_CP", &S_VV_S_CP);
  tree->SetBranchAddress("S_VV_P_CP", &S_VV_P_CP);
  tree->SetBranchAddress("S_VV_D_CP", &S_VV_D_CP);

  tree->SetBranchAddress("pdf_gen", &pdf_gen);

  for (std::size_t i = 0; i < n; ++i) {
    tree->GetEntry(i);

    Event event(*p4Kp, *p4pim, *p4Km, *p4pip, time, qtag, has_cache);

    if (Selection_masses(event.GetMassKpPim(), event.GetMassKmPip()) == false)
      continue;
    event.SetGenPDF(pdf_gen);

    event.SetEff(1.0);

    if (has_cache) {
      B_VV = {B_VV_S, B_VV_P, B_VV_D};
      B_VV_CP = {B_VV_S_CP, B_VV_P_CP, B_VV_D_CP};
      event.GetToSetVV().SetAmpsBarrier(B_VV, B_VV_CP);

      S_VV = {*S_VV_S, *S_VV_P, *S_VV_D};
      S_VV_CP = {*S_VV_S_CP, *S_VV_P_CP, *S_VV_D_CP};
      event.GetToSetVV().SetAmpsSpin(S_VV, S_VV_CP);
    }

    data.push_back(event);
  }

  file->Close();

  std::cout << data.size() << " toy events loaded" << std::endl;
}

void SaveData(const std::string &filename, const std::vector<Event> &data) {
  TFile *const file = new TFile(filename.c_str(), "RECREATE");

  TTree *const tree = new TTree("fitTree", "");

  // Declare branches
  double time;
  int qtag;

  TLorentzVector p4Kp, p4pim, p4Km, p4pip;

  double B_VV_S, B_VV_P, B_VV_D;
  double B_VV_S_CP, B_VV_P_CP, B_VV_D_CP;

  std::complex<double> S_VV_S, S_VV_P, S_VV_D;
  std::complex<double> S_VV_S_CP, S_VV_P_CP, S_VV_D_CP;

  double cos1, cos2, phi, m1, m2;

  double eff;

  double pdf_gen;

  // Set Branch adresses
  tree->Branch("time", &time);
  tree->Branch("qtag", &qtag);

  tree->Branch("p4Kp", &p4Kp);
  tree->Branch("p4pim", &p4pim);
  tree->Branch("p4Km", &p4Km);
  tree->Branch("p4pip", &p4pip);

  tree->Branch("eff", &eff);

  tree->Branch("B_VV_S", &B_VV_S);
  tree->Branch("B_VV_P", &B_VV_P);
  tree->Branch("B_VV_D", &B_VV_D);
  tree->Branch("B_VV_S_CP", &B_VV_S_CP);
  tree->Branch("B_VV_P_CP", &B_VV_P_CP);
  tree->Branch("B_VV_D_CP", &B_VV_D_CP);

  tree->Branch("S_VV_S", &S_VV_S);
  tree->Branch("S_VV_P", &S_VV_P);
  tree->Branch("S_VV_D", &S_VV_D);
  tree->Branch("S_VV_S_CP", &S_VV_S_CP);
  tree->Branch("S_VV_P_CP", &S_VV_P_CP);
  tree->Branch("S_VV_D_CP", &S_VV_D_CP);

  tree->Branch("CosTheta1", &cos1);
  tree->Branch("CosTheta2", &cos2);
  tree->Branch("Phi", &phi);
  tree->Branch("mKpPim", &m1);
  tree->Branch("mKmPip", &m2);

  tree->Branch("pdf_gen", &pdf_gen);

  for (unsigned int i = 0; i < data.size(); ++i) {
    const Event &event = data[i];

    time = event.GetTime();
    qtag = event.GetQtag();
    p4Kp = event.Getp4(CPConf::A)[0].Getp4ROOT();
    p4pim = event.Getp4(CPConf::A)[1].Getp4ROOT();
    p4Km = event.Getp4(CPConf::A)[2].Getp4ROOT();
    p4pip = event.Getp4(CPConf::A)[3].Getp4ROOT();

    B_VV_S = event.GetVV().GetAmpBarrier(SpinVV::S, CPConf::A);
    B_VV_P = event.GetVV().GetAmpBarrier(SpinVV::P, CPConf::A);
    B_VV_D = event.GetVV().GetAmpBarrier(SpinVV::D, CPConf::A);
    B_VV_S_CP = event.GetVV().GetAmpBarrier(SpinVV::S, CPConf::Abar);
    B_VV_P_CP = event.GetVV().GetAmpBarrier(SpinVV::P, CPConf::Abar);
    B_VV_D_CP = event.GetVV().GetAmpBarrier(SpinVV::D, CPConf::Abar);

    S_VV_S = event.GetVV().GetAmpSpin(SpinVV::S, CPConf::A);
    S_VV_P = event.GetVV().GetAmpSpin(SpinVV::P, CPConf::A);
    S_VV_D = event.GetVV().GetAmpSpin(SpinVV::D, CPConf::A);
    S_VV_S_CP = event.GetVV().GetAmpSpin(SpinVV::S, CPConf::Abar);
    S_VV_P_CP = event.GetVV().GetAmpSpin(SpinVV::P, CPConf::Abar);
    S_VV_D_CP = event.GetVV().GetAmpSpin(SpinVV::D, CPConf::Abar);

    m1 = (p4Kp + p4pim).M();
    m2 = (p4Km + p4pip).M();
    cos1 = event.GetCosTheta1();
    cos2 = event.GetCosTheta2();
    phi = event.GetChi();

    eff = event.GetEff();

    pdf_gen = event.GetGenPDF();

    tree->Fill();
  }

  tree->Write();
  file->Close();
}
