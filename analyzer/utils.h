#include <vector>

#include "TH1F.h"
#include "TLorentzVector.h"

using namespace std;

// deltaR:
double fdeltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaphi = fabs(phi1 - phi2);
  if (deltaphi > TMath::Pi())
    deltaphi = TMath::TwoPi() - deltaphi;
  double deltaeta = fabs(eta1 - eta2);
  double deltaR = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);
  return deltaR;
}

// deltaPhi:
double fdeltaPhi(double phi1, double phi2)
{
  double deltaphi = fabs(phi1 - phi2);
  if (deltaphi > TMath::Pi())
    deltaphi = TMath::TwoPi() - deltaphi;
  return deltaphi;
}

// Hemispheres:
vector<TLorentzVector> CombineJets(vector<TLorentzVector> myjets)
{
  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;
  //bool foundGood = false;
  int N_comb = 1;
  for(unsigned int i = 0; i < myjets.size(); i++){
    N_comb *= 2;
  }
  double M_min = 9999999999.0;
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
        j_temp1 += myjets[count];
      } else {
        j_temp2 += myjets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.M2()+j_temp2.M2();
    // smallest mass
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }  
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }
  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;
}


// MR
double CalcMR(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();
  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());
  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
  double mygamma = 1./sqrt(1.-mybeta*mybeta);
  //gamma times MRstar
  temp *= mygamma;
  return temp;
}


// MTR
double CalcMTR(TLorentzVector ja, TLorentzVector jb, TVector3 met){
  double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
  temp /= 2.;
  temp = sqrt(temp);
  return temp;
}

// MT
double CalcMT(TLorentzVector lepton, TLorentzVector pfmet){
  return sqrt( 2 * lepton.Pt() * pfmet.Pt() * ( 1 - cos( pfmet.Phi() - lepton.Phi() ) ) );
}

// Scalefactor for top pt reweighting
// Taken from this twiki: https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#Studies
// Values used are the ones for the 8 TeV, all combined measurement
double GetTopPtScaleFactor(double toppt){ // toppt is generator top(antitop) pt
  double a = 0.156;
  double b = -0.00137;
  return exp(a+b*toppt);
}

Double_t* getFixedBinEdges(int nbins, double mini, double maxi)
{
  Double_t* bin_edges = new Double_t[nbins+1];
  Double_t spacing = (maxi-mini)/nbins;
  for (int i = 0; i != nbins; ++i) {
    bin_edges[i] = mini + i*spacing;
    //cout << "bin edge " << i << " : " << bin_edges[i] << endl;
  }
  bin_edges[nbins] = maxi;

  return bin_edges;
}

Double_t* getVariableBinEdges(int num_entries, Double_t* tmp_array)
{
  Double_t* my_array = new Double_t[num_entries];
  for (int i = 0; i != num_entries; ++i) {
    my_array[i] = tmp_array[i];

    //cout << "bin edge " << i << " : " << my_array[i] << endl;
  }
   return my_array;
}

// Get efficiency from a 1D histogram
double geteff1D(TH1D* h, double x)
{
  double eff = 0.0;
  for (int i=1; i<h->GetNbinsX()+1; i++) {
    double xmin = h->GetXaxis()->GetBinLowEdge(i);
    double xmax = h->GetXaxis()->GetBinUpEdge(i);
    if (!(x >= xmin and x < xmax)) continue;
    eff = h->GetBinContent(i);
    break;
  }
  return eff;
}

// Get efficiency from a 2D histogram
double geteff2D(TH2D* h, double x, double y)
{
  double eff = 0.0;
  for (int i=1; i<h->GetNbinsX()+1; i++) {
    double xmin = h->GetXaxis()->GetBinLowEdge(i);
    double xmax = h->GetXaxis()->GetBinUpEdge(i);
    if (!(x >= xmin and x < xmax)) continue;
    for (int j=1; j<h->GetNbinsY()+1; j++) {
      double ymin = h->GetYaxis()->GetBinLowEdge(j);
      double ymax = h->GetYaxis()->GetBinUpEdge(j);
      if (y >= ymin && y < ymax) {
	eff = h->GetBinContent(i, j);
	break;
      }
    }
  }
  return eff;
}
