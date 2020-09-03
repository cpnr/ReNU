#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

void displayDetector()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // Read detector geometry file
  double rmax=0;
  std::map<int, std::array<double,3> > pmtIdToPosMap;
  std::ifstream fin("../data/pmtcoordinates_ID.dat");
  if ( true ) {
    int pmtId, dummy;
    double x, y, z;
    while ( !fin.eof() ) {
      fin >> pmtId;
      if ( fin.eof() ) break;
      fin >> x >> y >> z >> dummy;

      pmtIdToPosMap[pmtId] = {{x,y,z}};
      if ( z >= -2720 and z <= 2720 ) rmax = std::max(rmax, std::hypot(x, y));
    }
  }

  // PMT location
  TGraph* gUpper = new TGraph();
  TGraph* gLower = new TGraph();
  TGraph* gSide = new TGraph();

  double wmin = 1e9, wmax = -1e9, hmin = 1e9, hmax = -1e9;
  for ( auto iter = pmtIdToPosMap.begin(); iter != pmtIdToPosMap.end(); ++iter ) {
    const int pmtId = iter->first;
    const double x = iter->second[0];
    const double y = iter->second[1];
    const double z = iter->second[2];

    if ( z >= -2720 and z <= 2720 ) {
      const double r = std::hypot(x, y);
      const double phi = std::atan2(y, x);
      const double d = r*phi;
      gSide->SetPoint(gSide->GetN(), d, z);

      wmin = std::min(wmin, d);
      wmax = std::max(wmax, d);
      hmin = std::min(hmin, z);
      hmax = std::max(hmax, z);
    }
    else {
      if      ( z >=  2720 ) gUpper->SetPoint(gUpper->GetN(), x, y+z+rmax);
      else if ( z <= -2720 ) gLower->SetPoint(gLower->GetN(), x, y+z-rmax);

      wmin = std::min(wmin, x);
      wmax = std::max(wmax, x);
      hmin = std::min(hmin, y+z-rmax);
      hmax = std::max(hmax, y+z+rmax);
    }
  }
  wmax *= 1.1;
  wmin *= 1.1;
  hmax *= 1.1;
  hmin *= 1.1;

  TCanvas* c = new TCanvas("c", "c", 700, 700*(hmax-hmin)/(wmax-wmin));
  TH2D* hFrame = new TH2D("hFrame", "", 100, wmin, wmax, 100, hmin, hmax);
  hFrame->Draw();
  TGraph* grps[] = {gSide, gUpper, gLower};
  for ( auto g : grps ) {
    g->SetMarkerSize(1);
    g->SetMarkerStyle(kOpenCircle);
    g->SetEditable(false);
    g->Draw("Psame");
  }

}

