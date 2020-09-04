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
  DetHist();
  //~DetHist();
  void Fill(const int detId, const double w=1.0);
  void Draw(TString opt) {hDet_->Draw(opt);};
  int GetRegion(const double z) { return z < -2300 ? -1 : z > 2300 ? +1 : 0; }

private:
  std::unique_ptr<TH2Poly> hDet_;
  
};

DetHist::DetHist()
{
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

      for ( unsigned k=0; k<ys.size(); ++k ) ys[k] += 5000;
      hDet_->AddBin(xs.size(), &xs[0], &ys[0]);
      for ( unsigned k=0; k<ys.size(); ++k ) ys[k] -= 10000;
      hDet_->AddBin(xs.size(), &xs[0], &ys[0]);
    }
  }
}

void showLayoutFixed()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  DetHist* hDet = new DetHist;

  const double hmin = -8500, hmax = 8500, wmin = -8500, wmax = 8500;
  TCanvas* c = new TCanvas("c", "c", 700, 700*(hmax-hmin)/(wmax-wmin));
  TH2D* hFrame = new TH2D("hFrame", "", 100, wmin, wmax, 100, hmin, hmax);
  hFrame->Draw();
  hDet->Draw("same");
}

