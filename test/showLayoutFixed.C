#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <memory>

#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMath.h"

using namespace std;

class DetHist
{
public:
  DetHist(const std::string name);
  //~DetHist();
  void Fill(const int id, const double w=1.0);
  void SetContent(const int id, const double binContent);
  void Draw(TString opt) { hDet_->Draw(opt); }
  int GetRegionByZ(const double z) const { return z < -2300 ? -1 : z > 2300 ? +1 : 0; }
  int GetRegionById(const int pmtId) const { return GetRegionByZ(pmtIdToPosMap_.at(pmtId)[2]); }
  int GetN() const { return pmtIdToPosMap_.size(); }
  void SetRange(const double zmin, const double zmax) { hDet_->SetMinimum(zmin); hDet_->SetMaximum(zmax); }

private:
  std::unique_ptr<TH2Poly> hDet_;
  std::map<int, std::array<double,3> > pmtIdToPosMap_;
  const double yShiftOnCanvas_ = 5000;

};

DetHist::DetHist(const std::string name)
{
  // Read detector geometry file
  double rmax=0;
  std::ifstream fin("../data/pmtcoordinates_ID.dat");
  int pmtId, dummy;
  double x, y, z;
  while ( !fin.eof() ) {
    fin >> pmtId;
    if ( fin.eof() ) break;
    fin >> x >> y >> z >> dummy;

    pmtIdToPosMap_[pmtId] = {{x,y,z}};
  }

  // Build TH2Poly histogram
  const std::vector<double> dSide = {
     7627.27,  7017.09,  6406.91,  5796.73,  5186.55,  4576.36,  3966.18,
     3356   ,  2745.82,  2135.64,  1525.45,  915.273,  305.091, -305.091,
    -915.273, -1525.45, -2135.64, -2745.82, -3356   , -3966.18, -4576.36,
    -5186.55, -5796.73, -6406.91, -7017.09, -7627.27,
  };
  const std::vector<double> zSide = {
    2200, 1650, 1100, 550, 0, -550, -1100, -1650, -2200,
  };

  const std::vector<double> xSect = {
    550, 1062.52, 1062.52,
    1550.49, 1650, 1550.49,
    2032.53, 2181.18, 2181.18, 2032.53,
  };
  const std::vector<double> ySect = {
    0, 284.701, -284.701,
    564.333, 0, -564.333,
    841.904, 287.158, -287.158, -841.904,
  };

  // Shapes
  hDet_.reset(new TH2Poly());
  hDet_->SetName(name.c_str());

  // For side
  for ( unsigned i=0; i<dSide.size(); ++i ) {
    for ( unsigned j=0; j<zSide.size(); ++j ) {
      const double w = 220, h = 200;
      hDet_->AddBin(dSide[i]-w, zSide[j]-h, dSide[i]+w, zSide[j]+h);
    }
  }
  // For upper and lower caps
  for ( unsigned i=0; i<xSect.size(); ++i ) {
    const double x0 = xSect[i], y0 = ySect[i];
    const double w = 340, h = 250;
    for ( unsigned j=0; j<6; ++j ) {
      const double t = 2*TMath::Pi()*j/6;
      const double x =  x0*cos(t) + y0*sin(t);
      const double y = -x0*sin(t) + y0*cos(t);

      std::array<double,4> xs = {{x0-w, x0, x0+w, x0}};
      std::array<double,4> ys = {{y0, y0+h, y0, y0-h}};
      for ( unsigned k=0; k<xs.size(); ++k ) {
        const double newX =  xs[k]*cos(t) + ys[k]*sin(t);
        const double newY = -xs[k]*sin(t) + ys[k]*cos(t);
        xs[k] = newX;
        ys[k] = newY;
      }

      for ( unsigned k=0; k<ys.size(); ++k ) ys[k] += yShiftOnCanvas_;
      hDet_->AddBin(xs.size(), &xs[0], &ys[0]);
      for ( unsigned k=0; k<ys.size(); ++k ) ys[k] -= 2*yShiftOnCanvas_;
      hDet_->AddBin(xs.size(), &xs[0], &ys[0]);
    }
  }
}

void DetHist::Fill(const int id, const double w=1.0)
{
  const int region = GetRegionById(id);
  const auto pos = pmtIdToPosMap_[id];
  if ( region == 0 ) {
    const double r = std::hypot(pos[0], pos[1]);
    const double phi = std::atan2(pos[1], pos[0]);
    hDet_->Fill(r*phi, pos[2], w);
  }
  else {
    const double shift = yShiftOnCanvas_*GetRegionById(id);
    hDet_->Fill(pos[0], pos[1]+shift, w);
  }
}

void DetHist::SetContent(const int id, const double content)
{
  const int region = GetRegionById(id);
  const auto pos = pmtIdToPosMap_[id];
  if ( region == 0 ) {
    const double r = std::hypot(pos[0], pos[1]);
    const double phi = std::atan2(pos[1], pos[0]);
    const int b = hDet_->FindBin(r*phi, pos[2]);
    hDet_->SetBinContent(b, content);
  }
  else {
    const double shift = yShiftOnCanvas_*GetRegionById(id);
    const int b = hDet_->FindBin(pos[0], pos[1]+shift);
    hDet_->SetBinContent(b, content);
  }
}

void showLayoutFixed()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const double hmin = -8500, hmax = 8500, wmin = -8500, wmax = 8500;
  TCanvas* c = new TCanvas("c", "c", 2*300, 2*300*(hmax-hmin)/(wmax-wmin));
  c->Divide(2,2);
  //TH2D* hFrame = new TH2D("hFrame", "", 100, wmin, wmax, 100, hmin, hmax);
  //hFrame->Draw();

  DetHist* hDetTime1 = new DetHist("hDetTime1");
  DetHist* hDetTime2 = new DetHist("hDetTime2");
  DetHist* hDetCnt1 = new DetHist("hDetCnt1");
  DetHist* hDetCnt2 = new DetHist("hDetCnt2");
  //for ( int i=0; i<hDet->GetN(); ++i ) hDet->Fill(i, i+1);

  TFile* f = TFile::Open("../data/IBD_MC_for_ML.root");
  TTree* tree = (TTree*)f->Get("event");
  const int nEntries = tree->GetEntries();

  int b_nHits;
  const int max_b_nHits = 1000;
  int b_hitCounts[max_b_nHits], b_hitPMTIds[max_b_nHits];
  double b_hitTimes[max_b_nHits];
  tree->SetBranchAddress("photon_hits", &b_nHits);
  tree->SetBranchAddress("hit_count", &b_hitCounts);
  tree->SetBranchAddress("hit_time", &b_hitTimes);
  tree->SetBranchAddress("hit_pmt", &b_hitPMTIds);
  while ( true ) {
    int iEvent;
    cout << "Type event number to display ([0-" << (nEntries-1) << "]:";
    cin >> iEvent;
    if ( iEvent < 0 ) break;

    tree->GetEntry(iEvent);
    
    double t1Min = 1e9, t1Max = -1e9;
    double t2Min = 1e9, t2Max = -1e9;
    for ( int i=0; i<b_nHits; ++i ) {
      //cout << b_hitTimes[i] << ' ';
      const int id = b_hitPMTIds[i];
      const double t = b_hitTimes[i];
      const int cnt = b_hitCounts[i];
      if ( t < 1000 ) {
        hDetTime1->SetContent(id, t);
        t1Min = std::min(t1Min, t);
        t1Max = std::max(t1Max, t);
        hDetCnt1->Fill(id, cnt);
      }
      else {
        hDetTime2->SetContent(id, t);
        t2Min = std::min(t2Min, t);
        t2Max = std::max(t2Max, t);
        hDetCnt2->Fill(id, cnt);
      }
    }
    hDetTime1->SetRange(t1Min, t1Max);
    hDetTime2->SetRange(t2Min, t2Max);

    c->cd(1);
    hDetTime1->Draw("COLZ");
    c->cd(2);
    hDetTime2->Draw("COLZ");
    c->cd(3);
    hDetCnt1->Draw("COLZ");
    c->cd(4);
    hDetCnt2->Draw("COLZ");

    c->Update();
  }
}

