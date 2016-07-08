#include "Common.h"

//-- standard root includes                                                                                                                                                                                  
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TGraphErrors.h>

//-- standard c++ includes                                                                                                                                                                                   
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>

using namespace std;

Common::Common() { }

Common::~Common() { }

void Common::CheckHisto(TH1D* h) {
  
  Int_t NbinX = h->GetNbinsX();
  
  if(h->GetBinContent(0) != 0) cout<<"histo "<<h->GetName()<<" with title "<<h->GetTitle()<<" has underflow of "<<h->GetBinContent(0)<<" entries"<<endl;
  if(h->GetBinContent(NbinX+1) != 0) cout<<"histo "<<h->GetName()<<" with title "<<h->GetTitle()<<" has overflow of "<<h->GetBinContent(NbinX+1)<<" entries"<<endl;
  if(h->GetBinContent(0) == 0 && h->GetBinContent(NbinX+1) == 0) cout<<"histo "<<h->GetName()<<" with title "<<h->GetTitle()<<" has no underflow and no overflow"<<endl;

  return;
}

void Common::CheckHisto(TH2D* h) {
  
  Int_t NbinX = h->GetNbinsX();
  Int_t NbinY = h->GetNbinsY();
  
  double underflow_x = 0;
  double underflow_y = 0;
  
  double overflow_x = 0;
  double overflow_y = 0;
  
  for(int ibiny = 0; ibiny <= NbinY+1; ibiny++) {
    if(h->GetBinContent(0,ibiny) != 0) underflow_x+=h->GetBinContent(0,ibiny);
    if(h->GetBinContent(NbinX+1,ibiny) != 0) overflow_x+=h->GetBinContent(NbinX+1,ibiny);
  }
  
  for(int ibinx = 0; ibinx <= NbinX+1; ibinx++) {
    if(h->GetBinContent(ibinx,0) != 0)underflow_y+=h->GetBinContent(ibinx,0);
    if(h->GetBinContent(ibinx,NbinY+1) != 0) overflow_y+=h->GetBinContent(ibinx,NbinY+1);
  }
  
  if(underflow_x != 0) cout<<"histo "<<h->GetName()<<" with title "<<h->GetTitle()<<" has x underflow of "<<underflow_x<<endl;
  if(overflow_x != 0)  cout<<"histo "<<h->GetName()<<" with title "<<h->GetTitle()<<" has x overflow of "<<overflow_x<<endl;
  
  if(underflow_y != 0) cout<<"histo "<<h->GetName()<<" with title "<<h->GetTitle()<<" has y underflow of "<<underflow_y<<endl;
  if(overflow_y != 0)  cout<<"histo "<<h->GetName()<<" with title "<<h->GetTitle()<<" has y overflow of "<<overflow_y<<endl;
  
  return;
}

void::Common::CheckGraph(TGraphErrors* g) {
  cout<<"graph "<<g->GetName()<<" with title "<<g->GetTitle()<<" has been checked"<<endl;
}

void Common::SetAxisName(TH1D* h, TString xleg, TString yleg) {
  h->GetXaxis()->SetTitle(xleg);
  h->GetYaxis()->SetTitle(yleg);
}

TH1D* Common::MakeHisto(TString name, TString title, TString xleg, TString yleg, int nbin, double bmin, double bmax) {


  TH1D* h = new TH1D(name,title,nbin,bmin,bmax);

  h->GetXaxis()->SetTitle(xleg);
  h->GetYaxis()->SetTitle(yleg);

  h->SetMinimum(0);
  h->Sumw2();

  cout<<"create histo: "<<h->GetName()<<" with title: "<<h->GetTitle()<<endl;

  return h;
}

TH1D* Common::MakeHisto(TString name, TString title, TString xleg, TString yleg, int nbin, double* listbin) {

  TH1D* h = new TH1D(name,title,nbin,listbin);

  h->GetXaxis()->SetTitle(xleg);
  h->GetYaxis()->SetTitle(yleg);

  h->SetMinimum(0);
  h->Sumw2();

  cout<<"create histo: "<<h->GetName()<<" with title: "<<h->GetTitle()<<endl;

  return h;
}

TH2D* Common::MakeHisto(TString name, TString title, TString xlg, TString ylg, int nbinx, double bxmin, double bxmax, int nbiny, double bymin, double bymax) {


  TH2D* h = new TH2D(name,title,nbinx,bxmin,bxmax,nbiny,bymin,bymax);

  h->GetXaxis()->SetTitle(xlg);
  h->GetYaxis()->SetTitle(ylg);

  h->SetMinimum(0);
  h->Sumw2();

  cout<<"create histo: "<<h->GetName()<<" with title: "<<h->GetTitle()<<endl;

  return h;
}

TH2D* Common::MakeHisto(TString name, TString title, TString xleg, TString yleg, int nbinx, double* listbinx, int nbiny, double* listbiny) {

  TH2D* h = new TH2D(name,title,nbinx,listbinx,nbiny,listbiny);

  h->GetXaxis()->SetTitle(xleg);
  h->GetYaxis()->SetTitle(yleg);

  h->SetMinimum(0);
  h->Sumw2();

  cout<<"create histo: "<<h->GetName()<<" with title: "<<h->GetTitle()<<endl;

  return h;
}

void Common::MakeGraph(TGraphErrors* g, TString name, TString title, TString xleg, TString yleg) {

  g->SetName(name);
  g->SetTitle(title);

  g->GetXaxis()->SetTitle(xleg);
  g->GetYaxis()->SetTitle(yleg);
}

TProfile* Common::MakeProfile(TString name, TString title, TString xleg, TString yleg, int nbin, double bmin, double bmax) {

  TProfile* h = new TProfile(name,title,nbin,bmin,bmax);

  h->GetXaxis()->SetTitle(xleg);
  h->GetYaxis()->SetTitle(yleg);

  h->SetMinimum(0);
  h->Sumw2();

  cout<<"create profile: "<<h->GetName()<<" with title: "<<h->GetTitle()<<endl;

  return h;
}

void Common::ComputeRatio(double a,double b,double error_a,double error_b, double &ratio, double &error) {

  //-- ratio = a/b                                                                                                                                          
                                                
  ratio = 0;
  error = 0;

  if(b == 0) cout<<"you are going to divide by zero"<<endl;

  if(a*b != 0) {
    ratio = a/b;
    double error2 = TMath::Power(ratio,2)*(TMath::Power(error_a/a,2)+TMath::Power(error_b/b,2));
    error = TMath::Sqrt(error2);
  }

}

void Common::ComputeProduct(double a, double b, double error_a, double error_b, double &product, double &error) {

  product = 0;
  error = 0;

  product = a*b;
  
  if(a*b != 0) {
    double error2 = TMath::Power(product,2)*(TMath::Power(error_a/a,2)+TMath::Power(error_b/b,2));
    error = TMath::Sqrt(error2);
  }

}

void Common::ComputeHistoRatio(TH1D* hup, TH1D* hdown, TH1D& hratio) {

  for(int ibin = 1; ibin <= hup->GetNbinsX(); ibin++) {

    double up = hup->GetBinContent(ibin);
    double down = hdown->GetBinContent(ibin);

    double error_up = hup->GetBinError(ibin);
    double error_down = hdown->GetBinError(ibin);

    double ratio = 0;
    double error = 0;
    ComputeRatio(up,down,error_up,error_down,ratio,error);

    hratio.SetBinContent(ibin,ratio);
    hratio.SetBinError(ibin,error);
  }

}

void Common::ComputeHistoProduct(TH1D* ha, TH1D* hb, TH1D& hproduct) {

  for(int ibin = 1; ibin <= ha->GetNbinsX(); ibin++) {

    double a = ha->GetBinContent(ibin);
    double b = hb->GetBinContent(ibin);

    double error_a = ha->GetBinError(ibin);
    double error_b = hb->GetBinError(ibin);

    double product = 0;
    double error = 0;
    ComputeProduct(a,b,error_a,error_b,product,error);

    hproduct.SetBinContent(ibin,product);
    hproduct.SetBinError(ibin,error);
  }

}

void Common::ComputeOneMinusFake(TH1D* hFake, TH1D& hOneMinusFake) {

  for(int ibin = 1; ibin <= hFake->GetNbinsX(); ibin++) {

    double bin_content = 1-hFake->GetBinContent(ibin);
    double bin_error = hFake->GetBinError(ibin);

    hOneMinusFake.SetBinContent(ibin,bin_content);
    hOneMinusFake.SetBinError(ibin,bin_error);
  }
}

void Common::ComputeCFPt(TH1D* hOneMinusBkg, TH1D* hAcc, TH1D* hPurity, TH1D* hStability, TH1D& hCF) {

  for(int ibin = 1; ibin <= hOneMinusBkg->GetNbinsX(); ibin++) {

    double one_minus_bkg = hOneMinusBkg->GetBinContent(ibin);
    double acceptance = hAcc->GetBinContent(ibin);
    double purity = hPurity->GetBinContent(ibin);
    double stability = hStability->GetBinContent(ibin);

    double bkg_error = hOneMinusBkg->GetBinError(ibin);
    double acc_error = hAcc->GetBinError(ibin);
    double purity_error = hPurity->GetBinError(ibin);
    double stability_error = hStability->GetBinError(ibin);

    double CF = 0;
    double CF_error2 = 0;
    double CF_error = 0;

    if(acceptance*stability != 0) 
      CF = (one_minus_bkg*purity)/(acceptance*stability);
    
    if(one_minus_bkg*acceptance*purity*stability != 0) 
      CF_error2 = CF*CF*(TMath::Power(bkg_error/one_minus_bkg,2) + TMath::Power(acc_error/acceptance,2) 
			+ TMath::Power(purity_error/purity,2) + TMath::Power(stability_error/stability,2));
    
    CF_error = TMath::Sqrt(CF_error2);
  
    hCF.SetBinContent(ibin,CF);
    hCF.SetBinError(ibin,CF_error);
  }
}


int Common::GetPhiMapping(int sector) {

  int ibin = 0;

  if(sector == 1) ibin = 13;
  if(sector == 2) ibin = 14;
  if(sector == 3) ibin = 15;
  if(sector == 4) ibin = 16;

  if(sector == 5) ibin = 1;
  if(sector == 6) ibin = 2;
  if(sector == 7) ibin = 3;
  if(sector == 8) ibin = 4;
  if(sector == 9) ibin = 5;
  if(sector == 10) ibin = 6;
  if(sector == 11) ibin = 7;
  if(sector == 12) ibin = 8;
  if(sector == 13) ibin = 9;
  if(sector == 14) ibin = 10;
  if(sector == 15) ibin = 11;
  if(sector == 16) ibin = 12;

  return(ibin);
}

void Common::GetHistoContent(TH1D* histo, double &content, double &error) {

  double error2 = 0;

  content = 0;
  error = 0;

  for(int ibin = 1; ibin <= histo->GetNbinsX(); ibin++) {
    content+=histo->GetBinContent(ibin);
    error2+=TMath::Power(histo->GetBinError(ibin),2);
  }

  error = TMath::Sqrt(error2);
}

double Common::GetDeltaR(double phi_trackjet, double eta_trackjet, double phi_genjet, double eta_genjet) {

  double DeltaEta = eta_trackjet - eta_genjet;
  double DeltaPhi = phi_trackjet - phi_genjet;

  double DeltaR2 = TMath::Power(DeltaEta,2) + TMath::Power(DeltaPhi,2);
  double DeltaR = TMath::Sqrt(DeltaR2);

  return(DeltaR);
}

int Common::GetSectorSubtracted(int isector_max, int shift) {

  //-- shift from -1 to -7
  int result = 0;

  if(isector_max - shift >= 0) 
    result = isector_max - shift;

  if(isector_max - shift < 0)
    result = isector_max - shift + 16;

  return(result);
}


int Common::GetSectorAdded(int isector_max, int shift) {

  //-- shift from +1 to +8
  int result = 0;

  if(isector_max + shift <= 15)
    result = isector_max + shift ;


  if(isector_max + shift > 15)
    result = isector_max + shift - 16;

  return(result);
}


double Common::GetDeltaPhi(double phi1, double phi2) {
 
  double result = phi1 - phi2;
  
  if(result > M_PI) result -= 2*M_PI;
  
  if(result <= -M_PI) result += 2*M_PI;
  
  return(result);
}

double Common::Rad2Deg(double phi) {
  //-- cout<<"phi = "<<phi<<" radians - ";
  double rad2deg = 180/M_PI;
  if(phi < 0) phi+=2*M_PI;
  phi*=rad2deg;
  //-- cout<<"phi = "<<phi<<" degrees"<<endl;
  return(phi);
}
  
TF1* Common::FitGauss(TH1D* histo,double min, double max, double DeltaM) {
  
  TF1* f = new TF1("fGauss","gaus",min,max);
  f->SetParNames("N","mean","sigma");
  f->SetLineColor(2);
  f->SetLineWidth(2);

  Int_t ibinmax = histo->GetMaximumBin();
  Float_t xmax = histo->GetBinLowEdge(ibinmax) + 0.5*histo->GetBinWidth(ibinmax);
  Float_t xlowedge = xmax - DeltaM*histo->GetRMS();
  Float_t xupedge = xmax + DeltaM*histo->GetRMS();
  f->SetRange(xlowedge,xupedge);
  histo->Fit("fGauss","RQNO");

  xmax = f->GetParameter(1);
  xlowedge = xmax - DeltaM*f->GetParameter(2);
  xupedge = xmax + DeltaM*f->GetParameter(2);
  f->SetRange(xlowedge,xupedge);
  histo->Fit("fGauss","RQNO");
  return f;
}

bool Common::IsGoodTrack(int quality,double pt, double pterror, double dxy, double dxyerror, double dz, double dzerror, int npixelhit) {

  bool goodtrack = true;

  //-- high quality
  if(quality != 2) goodtrack = false; 
  
  //-- pterror_pt <= 10
  double pterror_pt = 1000;
  if(pt != 0) pterror_pt = 100*pterror/pt;
  if(pterror_pt > 10) goodtrack = false;

  //-- dxysigmaxy <= 3
  double dxysigmaxy = 1000;
  if(dxyerror != 0) dxysigmaxy = dxy/dxyerror;
  if(dxysigmaxy > 3) goodtrack = false;

  //-- dzsigmaz <= 3
  double dzsigmaz = 1000;  
  if(dzerror != 0) dzsigmaz = dz/dzerror;
  if(dzsigmaz > 3) goodtrack = false;

  //-- n pixel hits >= 3
  if(npixelhit < 3) goodtrack = false;

  return(goodtrack);
}

double Common::GetVertexWeight(double Vz,bool* AllMC, int iteration) {

  //-- for(int i = 0; i < 3; i++) cout<<"all MC "<<i<<": "<<AllMC[i]<<endl;
  //-- getchar();

  double weight = 1;
  
  //-- normalisation
  double norma = 1;

  double norma_pythia = 1.00288;
  double norma_herwig = 1.00290;
  double norma_epos = 1.00259;
  
  if(AllMC[0] == true) norma = norma_pythia; //-- PYTHIA8
  if(AllMC[1] == true) norma = norma_herwig; //-- HERWIGPP
  if(AllMC[2] == true) norma = norma_epos;   //-- EPOS

  //-- vz range
  double Vzmin = -20;
  double Vzmax = 20;

  //-- monte carlo factor
  double factor_MC = 1;

  TF1* fweight_MC = new TF1("fweight_MC","gaus",Vzmin,Vzmax);
  double para_MC[3];

  //-- PYTHIA8
  if(AllMC[0] == true) {
    para_MC[0] = 12262.1;
    para_MC[1] = -1.08942;
    para_MC[2] = 5.2694;
  }

  //-- PYTHIA8
  //-- constant = 12262.1 +/- 15.7418
  //-- mean = -1.08942 +/- 0.00552802
  //-- sigma = 5.2694 +/- 0.00392203
  
  //-- HERWIGPP
  if(AllMC[1] == true) {
    para_MC[0] = 12310.7;
    para_MC[1] = -1.0712;
    para_MC[2] = 5.24863;
  }

  //-- HERWIGPP
  //-- constant = 12310.7 +/- 16.1184
  //-- mean = -1.0712 +/- 0.00561567
  //-- sigma = 5.24863 +/- 0.00398296

  //-- EPOS
  if(AllMC[2] == true) {
    para_MC[0] = 12268.8;
    para_MC[1] = -1.0829;
    para_MC[2] = 5.26641;
  }

  //-- EPOS
  //-- constant = 12268.8 +/- 19.4202
  //-- mean = -1.0829 +/- 0.00679314
  //-- sigma = 5.26641 +/- 0.00485983

  fweight_MC->SetParameters(para_MC);
  factor_MC = fweight_MC->Eval(Vz);

  //-- data factor
  double factor_data = 1;

  TF1* fweight_data = new TF1("fweight_data","gaus",Vzmin,Vzmax);
  double para_data[3];

  para_data[0] = 15337.5;
  para_data[1] = -1.67589;
  para_data[2] = 4.19475;

  //-- data
  //-- constant = 15337.5 +/- 28.2037
  //-- mean = -1.67589 +/- 0.00661041
  //-- sigma = 4.19475 +/- 0.00398527

  fweight_data->SetParameters(para_data);
  factor_data = fweight_data->Eval(Vz);

  //-- final weight
  if(iteration == 0) weight = 1;
  if(iteration == 1) weight = factor_data/factor_MC;
  if(iteration == 2) weight = norma*factor_data/factor_MC;
  if(!(Vz > Vzmin && Vz < Vzmax)) weight = 1;

  delete fweight_MC;
  delete fweight_data;

  return(weight);
}

void Common::GetPUWeight(bool* AllMC,double* weight,int iteration) {

  //-- iteration == 4 is the last one
  
  //-- PYTHIA8 AllMC[0]
  //-- HERWIGPP AllMC[1]
  //-- EPOS AllMC[2]     
  
  int nbin = 10;
  
  for(int ibin = 0; ibin < nbin; ibin++) weight[ibin] = 1;
  
  std::vector<double> weightPU_temp;  
  std::vector< std::vector<double> > weightPU;  
  
  //-- iteration 0
  for(int ibin = 0; ibin < nbin; ibin++) weightPU_temp.push_back(1);
  weightPU.push_back(weightPU_temp);
  weightPU_temp.clear();              
  
  //-- next iteration
  for(int iter = 0; iter < iteration && iteration > 0; iter++) {
    //-- retrieve file
    TString file_name = "Correction-PU/Correction-PU-";
    file_name+=TString::Format("%d",iter+1);
    file_name+=".root";
    TFile file_Correction(file_name);
    
    cout<<"file with PU correction will be open: "<<file_name<<endl;
    
    //-- retrieve histo
    TH1D *hCorrection = NULL;
    
    for(int i = 0; i < 3; i++) {
      TString histo_name = "histo_weight_";
      histo_name+=TString::Format("%d",i+1);
      if(AllMC[i] == true) hCorrection = (TH1D*) file_Correction.Get(histo_name); 
    }
    
    //-- retrieve weights
    for(int ibin = 0; ibin < nbin; ibin++) weightPU_temp.push_back(hCorrection->GetBinContent(ibin+1));
    weightPU.push_back(weightPU_temp);  
    weightPU_temp.clear();              
  }
  
  cout<<endl<<"iteration = "<<iteration<<endl<<endl;
  
  for(int iter = 0; iter < iteration+1; iter++) {
    for(int ibin = 0; ibin < nbin; ibin++) {
      cout<<"pu weight "<<iter<<" - "<<ibin+1<<" vertex = "<<weightPU.at(iter).at(ibin)<<endl;
      weight[ibin]*=weightPU.at(iter).at(ibin);
    }
    cout<<endl;
  }
  
  for(int ibin = 0; ibin < nbin; ibin++) 
    cout<<"pu weight global - "<<ibin+1<<" vertex = "<<weight[ibin]<<endl;
  
  //-- normalisation
  if(iteration == 4) {
    double norma = 1;
    
    double norma_pythia = 0.961774;
    double norma_herwig = 0.96349;
    double norma_epos = 0.962311;
    
    if(AllMC[0] == true) norma = norma_pythia; //-- PYTHIA8                            
    if(AllMC[1] == true) norma = norma_herwig; //-- HERWIGPP                                                          
    if(AllMC[2] == true) norma = norma_epos;   //-- EPOS    
    
    for(int ibin = 0; ibin < nbin; ibin++) weight[ibin]*=norma;
  }
  
}
   
bool Common::isProton(bool* AllMC, int id) {

  bool isproton = false;

  //-- PYTHIA8
  if(AllMC[0] && TMath::Abs(id) == 2212) isproton = true;

  //-- HERWIGPP
  if(AllMC[1] && TMath::Abs(id) == 2212) isproton = true;

  //-- EPOS
  if(AllMC[2] && TMath::Abs(id) == 2212) isproton = true;

  return(isproton);
}

bool Common::isPion(bool* AllMC, int id) {

  bool ispion = false;

  //-- PYTHIA8						 
  if(AllMC[0] && TMath::Abs(id) == 211) ispion = true;  
							     
  //-- HERWIGPP						 
  if(AllMC[1] && TMath::Abs(id) == 211) ispion = true;
							 
  //-- EPOS						 
  if(AllMC[2] && TMath::Abs(id) == 211) ispion = true;

  return(ispion);                                        
}

bool Common::isKaon(bool* AllMC, int id) {

  bool iskaon = false;					 
							 
  //-- PYTHIA8						 
  if(AllMC[0] && TMath::Abs(id) == 321) iskaon = true;   
							 
  //-- HERWIGPP						 
  if(AllMC[1] && TMath::Abs(id) == 321) iskaon = true;	 
							 
  //-- EPOS						 
  if(AllMC[2] && TMath::Abs(id) == 321) iskaon = true;	 
							 
  return(iskaon);                                        
}

