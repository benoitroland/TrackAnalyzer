#include "TrackAnalyzer.h"

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
#include <string>
#include <cstring>

#define debugMCh 0
#define HFmatching 0
#define HFleading 0
#define BScheck 0
#define GoodVertex 0
#define PUinfo 0
#define PUdebug 0
#define VertexDebug 0
#define HitDebug 0
#define MatchingDebug 0

TrackAnalyzer::TrackAnalyzer() { }

TrackAnalyzer::~TrackAnalyzer() { }

void TrackAnalyzer::Loop(TString inputdir,TObjArray* filelist,TString type, TString selection, TString file_name) {

  using namespace std;

  const double MY_PI = TMath::Pi();

  //-- type 
  bool isData = false;				       
  bool isPYTHIA8 = false;
  bool isHERWIGPP = false;
  bool isEPOS = false;

  if(strcmp(type,"data-0T") == 0) isData = true;	       
  if(strcmp(type,"data-38T") == 0) isData = true;
  if(strcmp(type,"data-HFRereco-38T") == 0) isData = true;
  if(strcmp(type,"data-PromptReco-38T") == 0) isData = true;

  if(strcmp(type,"PYTHIA8-CUETP8M1-38T-PU1p3") == 0) isPYTHIA8 = true;
  if(strcmp(type,"HERWIGPP-CUETHS1-38T-PU1p3") == 0) isHERWIGPP = true;
  if(strcmp(type,"EPOS-38T-PU1p3") == 0) isEPOS = true;  

  bool AllMC[3];
  for(int i = 0; i < 3; i++) AllMC[i] = false;
  if(isPYTHIA8 == true) AllMC[0] = true;
  if(isHERWIGPP == true) AllMC[1] = true;
  if(isEPOS == true) AllMC[2] = true;

  bool pu_reweight = false;
  bool pixel_cut = false;
  bool vz_reweight = false; 
  bool combined_cut = false;

  double HFscaling = 0;
				 
  if(strcmp(selection,"no_pu_reweight") == 0) {  
   pu_reweight = false; 		       
   pixel_cut = true;   		       
   vz_reweight = true; 		       

   combined_cut = false;

   if(isData) HFscaling = 1;
   if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }                                            

  if(strcmp(selection,"no_pixel_cut") == 0) {
    pu_reweight = true; 
    pixel_cut = false;   
    vz_reweight = true; 

    combined_cut = false;

    if(isData) HFscaling = 1;				
    if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }

  if(strcmp(selection,"no_vz_reweight") == 0) {
   pu_reweight = true; 
   pixel_cut = true;   
   vz_reweight = false; 
   
   combined_cut = false;

   if(isData) HFscaling = 1;				   
   if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }

  if(strcmp(selection,"nominal") == 0) {
    pu_reweight = true;
    pixel_cut = true;
    vz_reweight = true;

    combined_cut = false;

   if(isData) HFscaling = 1;				
   if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }                                               

  if(strcmp(selection,"HFup") == 0) {
    pu_reweight = true;
    pixel_cut = true;
    vz_reweight = true;

    combined_cut = false;

    if(isData) HFscaling = 1.20;				
    if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }
 
  if(strcmp(selection,"HFdown")  == 0) {
    pu_reweight = true;
    pixel_cut = true;
    vz_reweight = true;

    combined_cut = false;

    if(isData) HFscaling = 0.80;			
    if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }

  if(strcmp(selection,"combined") == 0) {
    pu_reweight = true;
    pixel_cut = false;
    vz_reweight = true;

    combined_cut = true;

   if(isData) HFscaling = 1;				
   if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }                                               

  cout<<endl;
  cout<<"type = "<<type<<endl;
  cout<<"selection = "<<selection<<endl<<endl;
  cout<<"pu reweight = "<<pu_reweight<<endl;
  cout<<"pixel cut = "<<pixel_cut<<endl;
  cout<<"combined cut = "<<combined_cut<<endl;
  cout<<"vz reweight = "<<vz_reweight<<endl<<endl;
  if(isData) cout<<"data - HF scaling = "<<HFscaling<<endl;	  
  if(isPYTHIA8) cout<<"PYTHIA8 - HF scaling = "<<HFscaling<<endl;	  
  if(isHERWIGPP) cout<<"HERWIGPP - HF scaling = "<<HFscaling<<endl;	  
  if(isEPOS) cout<<"EPOS - HF scaling = "<<HFscaling<<endl;	  
  cout<<endl;
  cout<<"press enter to continue"<<endl;
  getchar();
  
  //-- event
  long double Nevt_tot = 0;
  long double Nevt_tot_no_weight = 0;     
  long double Nevt_tot_no_pu_weight = 0;  
  long double Nevt_tot_no_vz_weight = 0;

  long double Ndata_sel = 0;
  long double Ndata_sel_no_weight = 0;    
  long double Ndata_sel_no_pu_weight = 0; 
  long double Ndata_sel_no_vz_weight = 0;

  long double Nmoca_reco_sel = 0;
  long double Nmoca_reco_sel_no_weight = 0;     
  long double Nmoca_reco_sel_no_pu_weight = 0;  
  long double Nmoca_reco_sel_no_vz_weight = 0;

  //double Nmoca_gen_sel = 0;
  
  double Nmoca_reco_sel_HFplus10 = 0;
  double Nmoca_reco_sel_HFplus20 = 0;

  double Nmoca_reco_sel_HFminus10 = 0;
  double Nmoca_reco_sel_HFminus20 = 0;

  //-- histo label and title

  TString label,title;
  TString xleg, yleg;

  //-- generated level information (no pT cut, no reweighting)											       

  TH1D *hChPt = MakeHisto("hChPt","charged particle pT","pT","N charged particles",300,0,30);
  TH1D *hChEta = MakeHisto("hChEta","charged particle eta","eta","N charged particles",200,-2.5,2.5);
  TH1D *hChPhi = MakeHisto("hChPhi","charged particle phi","phi","N charged particles",200,-TMath::Pi(),TMath::Pi());
  TH1D *hChNum = MakeHisto("hChNum","charged particle multiplicity","multiplicity","N events",500,0,500);                                    

  TH1D *hChLeadingPt = MakeHisto("hChLeadingPt","leading charged particle pT","leading pT","N events",300,0,30);			     
  TH1D *hChLeadingEta = MakeHisto("hChLeadingEta","leading charged particle eta","leading eta","N events",200,-2.5,2.5);		     
  TH1D *hChLeadingPhi = MakeHisto("hChLeadingPhi","leading charged particle phi","leading phi","N events",200,-TMath::Pi(),TMath::Pi());

  TH1D *h500ChPt = MakeHisto("h500ChPt","charged particle pT","pT","N charged particles",300,0,30);						
  TH1D *h500ChEta = MakeHisto("h500ChEta","charged particle eta","eta","N charged particles",200,-2.5,2.5);					
  TH1D *h500ChPhi = MakeHisto("h500ChPhi","charged particle phi","phi","N charged particles",200,-TMath::Pi(),TMath::Pi());			
  TH1D *h500ChNum = MakeHisto("h500ChNum","charged particle multiplicity","multiplicity","N events",500,0,500);                               
                                                                                                                                      
  TH1D *h500ChLeadingPt = MakeHisto("h500ChLeadingPt","leading charged particle pT","leading pT","N events",300,0,30);			
  TH1D *h500ChLeadingEta = MakeHisto("h500ChLeadingEta","leading charged particle eta","leading eta","N events",200,-2.5,2.5);		
  TH1D *h500ChLeadingPhi = MakeHisto("h500ChLeadingPhi","leading charged particle phi","leading phi","N events",200,-TMath::Pi(),TMath::Pi());

  TH1D *hHadPt = MakeHisto("hHadPt","charged hadron pT","pT","N charged hadrons",300,0,30);						     
  TH1D *hHadEta = MakeHisto("hHadEta","charged hadron eta","eta","N charged hadrons",200,-2.5,2.5);					     
  TH1D *hHadPhi = MakeHisto("hHadPhi","charged hadron phi","phi","N charged hadrons",200,-TMath::Pi(),TMath::Pi());			     
  TH1D *hHadNum = MakeHisto("hHadNum","charged hadron multiplicity","multiplicity","N events",500,0,500);                                    

  //-- matching vertex track - all track

  TH1D *hCheckDeltaR = MakeHisto("hCheckDeltaR","Check Delta R","Delta R","N tracks",500,0,0.05);
  TH1D *hCheckDeltapT = MakeHisto("hCheckDeltapT","Check Delta pT","Delta pT","N tracks",2000,-10,10);
  TH1D *hCheckDeltaNum = MakeHisto("hCheckDeltaNum","Check Delta Num","Delta Num","N events",11,-5.5,5.5);

  //-- selection
  TH1D *hselection_data = MakeHisto("hselection_data","Event Selection Data","step","N evts",6,0.5,6.5);
  TH1D *hselection_moca_reco = MakeHisto("hselection_moca_reco","Event Selection MC reco","step","N evts",3,0.5,3.5);
  cout<<endl;

  //-- HF energy cut

  const int nHFcut = 3;		    
  double EHFcut[nHFcut] = {0.5,3,5};

  //-- HF tower segmentation
  //-- http://home.fnal.gov/~ratnikov/lpc_jetmet/calorimetry/tower_type.txt
  
  const int nHFtower = 13;
  const int nHFtower_10 = 11;
  const int nHFtower_20 = 2;

  double HFplus_eta_min = 2.853;
  double HFplus_eta_max = 5.191;
  double HFminus_eta_min = -5.191;
  double HFminus_eta_max = -2.853;

  double HF_eta_width[nHFtower] = {0.111,0.175,0.175,0.175,0.175,0.175,0.174,0.178,0.172,0.175,0.178,0.173,0.302};
  double HF_phi_width[nHFtower] = {10,10,10,10,10,10,10,10,10,10,10,20,20};

  //-- HF tower plus - eta segmentation - all values
  double value_HFplus_eta_edge[nHFtower+1];
  double value_HFplus_eta_edge_10[nHFtower_10+1];
  double value_HFplus_eta_edge_20[nHFtower_20+1];

  value_HFplus_eta_edge[0] = HFplus_eta_min;  
  for(int ieta = 0; ieta < nHFtower; ieta++)
    value_HFplus_eta_edge[ieta+1] = value_HFplus_eta_edge[ieta] + HF_eta_width[ieta];  

  value_HFplus_eta_edge_10[0] = HFplus_eta_min;  
  for(int ieta = 0; ieta < nHFtower_10; ieta++)
    value_HFplus_eta_edge_10[ieta+1] = value_HFplus_eta_edge_10[ieta] + HF_eta_width[ieta];  

  value_HFplus_eta_edge_20[0] =  value_HFplus_eta_edge[nHFtower-2];
  for(int ieta = 0; ieta < nHFtower_20; ieta++)
    value_HFplus_eta_edge_20[ieta+1] = value_HFplus_eta_edge_20[ieta] + HF_eta_width[nHFtower-2+ieta]; 
  
  //-- HF tower plus - eta segmentation - (min,max) for each bin
  std::vector< std::vector<double> > HFplus_eta_edge;
 
  for(int ieta = 0; ieta < nHFtower; ieta++) {
    std::vector<double> temp_HFplus_eta_edge;
    temp_HFplus_eta_edge.push_back(value_HFplus_eta_edge[ieta]);
    temp_HFplus_eta_edge.push_back(value_HFplus_eta_edge[ieta+1]);
    HFplus_eta_edge.push_back(temp_HFplus_eta_edge);
  }
 
  //-- HF tower plus - phi segmentation - all values
  std::vector< std::vector<double> > HFplus_phi_edge;
  
  for(int ieta = 0; ieta < nHFtower; ieta++){
    int nphi = 360/HF_phi_width[ieta] + 1;

    std::vector<double> temp_HFplus_phi_edge;
    for(int iphi = 0; iphi < nphi; iphi++) temp_HFplus_phi_edge.push_back(iphi*HF_phi_width[ieta]);

    HFplus_phi_edge.push_back(temp_HFplus_phi_edge);
  }
  
  const int nphi_10 = 36;
  double value_HFplus_phi_edge_10[nphi_10+1];
  for(int iphi = 0; iphi < nphi_10+1; iphi++) value_HFplus_phi_edge_10[iphi] = iphi*10;
  
  const int nphi_20 = 18;
  double value_HFplus_phi_edge_20[nphi_20+1];
  for(int iphi = 0; iphi < nphi_20+1; iphi++) value_HFplus_phi_edge_20[iphi] = iphi*20;

  //-- HF tower minus - eta segmentation - all values
  double value_HFminus_eta_edge[nHFtower+1];
  double value_HFminus_eta_edge_10[nHFtower_10+1];
  double value_HFminus_eta_edge_20[nHFtower_20+1];

  for(int ieta = 0; ieta < nHFtower+1; ieta++)
    value_HFminus_eta_edge[ieta] = - value_HFplus_eta_edge[nHFtower-ieta];

  for(int ieta = 0; ieta < nHFtower_10+1; ieta++)
    value_HFminus_eta_edge_10[ieta] = - value_HFplus_eta_edge_10[nHFtower_10-ieta];

  for(int ieta = 0; ieta < nHFtower_20+1; ieta++)
    value_HFminus_eta_edge_20[ieta] = - value_HFplus_eta_edge_20[nHFtower_20-ieta];

  //-- HF tower minus - eta segmentation - (min,max) for each bin
  std::vector< std::vector<double> > HFminus_eta_edge;
  
  for(unsigned int ieta = 0; ieta < HFplus_eta_edge.size(); ieta++){
    std::vector<double> temp_HFminus_eta_edge;

    temp_HFminus_eta_edge.push_back(-HFplus_eta_edge.at(nHFtower-1-ieta).at(1));
    temp_HFminus_eta_edge.push_back(-HFplus_eta_edge.at(nHFtower-1-ieta).at(0));

    HFminus_eta_edge.push_back(temp_HFminus_eta_edge);
  }

  //-- HF tower minus - phi segmentation - all values
  std::vector< std::vector<double> > HFminus_phi_edge;
  
  for(int ieta = 0; ieta < nHFtower; ieta++){
    std::vector<double> temp_HFminus_phi_edge;

    int nphi = HFplus_phi_edge.at(nHFtower-1-ieta).size();
    for(int iphi = 0; iphi < nphi; iphi++) temp_HFminus_phi_edge.push_back(HFplus_phi_edge.at(nHFtower-1-ieta).at(iphi));

    HFminus_phi_edge.push_back(temp_HFminus_phi_edge);
  }

  double value_HFminus_phi_edge_10[nphi_10+1];
  for(int iphi = 0; iphi < nphi_10+1; iphi++) value_HFminus_phi_edge_10[iphi] = iphi*10;
  
  double value_HFminus_phi_edge_20[nphi_20+1];
  for(int iphi = 0; iphi < nphi_20+1; iphi++) value_HFminus_phi_edge_20[iphi] = iphi*20;

  //-- HF tower plus - check
  
  cout<<endl;
  cout<<"    -- HF tower plus segmentation -- "<<endl<<endl;
  
  for(int ieta = 0; ieta < nHFtower+1; ieta++) cout<<value_HFplus_eta_edge[ieta]<<" ";
  cout<<endl<<endl;

  for(int ieta = 0; ieta < nHFtower_10+1; ieta++) cout<<value_HFplus_eta_edge_10[ieta]<<" ";
  cout<<"have 10 degrees phi segmentation"<<endl;                                         
  for(int iphi = 0; iphi < nphi_10+1; iphi++) cout<<value_HFplus_phi_edge_10[iphi]<<" ";
  cout<<endl<<endl;

  for(int ieta = 0; ieta < nHFtower_20+1; ieta++) cout<<value_HFplus_eta_edge_20[ieta]<<" "; 
  cout<<"have 20 degrees phi segmentation"<<endl;
  for(int iphi = 0; iphi < nphi_20+1; iphi++) cout<<value_HFplus_phi_edge_20[iphi]<<" ";
  cout<<endl<<endl;
  
  for(int ieta = 0; ieta < nHFtower; ieta++) {
    cout<<"HF plus tower "<<ieta+1<<": "<<HFplus_eta_edge.at(ieta).at(0)<<" < eta < "<<HFplus_eta_edge.at(ieta).at(1)<<endl;
  
    cout<<"phi segmentation: ";
    int nphi = HFplus_phi_edge.at(ieta).size();
    
    for(int iphi = 0; iphi < nphi; iphi++) {
      cout<<HFplus_phi_edge.at(ieta).at(iphi)<<" ";
    }
    
    cout<<endl<<endl;
  }

  //-- getchar();

  //-- HF tower minus - check

  cout<<"    -- HF tower minus segmentation -- "<<endl<<endl;

  for(int ieta = 0; ieta < nHFtower+1; ieta++) cout<<value_HFminus_eta_edge[ieta]<<" ";
  cout<<endl<<endl;
  
  for(int ieta = 0; ieta < nHFtower_20+1; ieta++) cout<<value_HFminus_eta_edge_20[ieta]<<" "; 
  cout<<"have 20 degrees phi segmentation"<<endl;
  for(int iphi = 0; iphi < nphi_20+1; iphi++) cout<<value_HFminus_phi_edge_20[iphi]<<" ";
  cout<<endl<<endl;

  for(int ieta = 0; ieta < nHFtower_10+1; ieta++) cout<<value_HFminus_eta_edge_10[ieta]<<" ";
  cout<<"have 10 degrees phi segmentation"<<endl;                                         
  for(int iphi = 0; iphi < nphi_10+1; iphi++) cout<<value_HFminus_phi_edge_10[iphi]<<" ";
  cout<<endl<<endl;

  for(int ieta = 0; ieta < nHFtower; ieta++) {
    cout<<"HF minus tower "<<ieta+1<<": "<<HFminus_eta_edge.at(ieta).at(0)<<" < eta < "<<HFminus_eta_edge.at(ieta).at(1)<<endl;
    
    cout<<"phi segmentation: ";
    int nphi = HFminus_phi_edge.at(ieta).size();
    
    for(int iphi = 0; iphi < nphi; iphi++) {
      cout<<HFminus_phi_edge.at(ieta).at(iphi)<<" ";
    }

    cout<<endl<<endl;
  }

  cout<<endl;
  //-- getchar();
  
  //-- HF tower occupancy
  const int nHFOccupancyCut = 5;		    
  double EHFOccupancyCut[nHFOccupancyCut] = {0.5,3,5,10,20};

  //-- HF plus tower occupancy
  TH2D* hHFplus10Occupancy[nHFOccupancyCut];
  TH2D* hHFplus20Occupancy[nHFOccupancyCut];

  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    label = "hHFplus10Occupancy_" + TString::Format("%d",ihf+1);
    title = "HF plus 10 degrees tower occupancy for E particle > " + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "#eta";
    yleg = "#varphi";
    hHFplus10Occupancy[ihf] = MakeHisto(label,title,xleg,yleg,nHFtower_10,value_HFplus_eta_edge_10,nphi_10,value_HFplus_phi_edge_10);

    label = "hHFplus20Occupancy_" + TString::Format("%d",ihf+1);
    title = "HF plus 20 degrees tower occupancy for E particle > " + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "#eta";
    yleg = "#varphi";
    hHFplus20Occupancy[ihf] = MakeHisto(label,title,xleg,yleg,nHFtower_20,value_HFplus_eta_edge_20,nphi_20,value_HFplus_phi_edge_20);
  }

  cout<<endl<<endl;

  //-- HF minus tower occupancy
  TH2D* hHFminus10Occupancy[nHFOccupancyCut];
  TH2D* hHFminus20Occupancy[nHFOccupancyCut];

  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    label = "hHFminus10Occupancy_" + TString::Format("%d",ihf+1);
    title = "HF minus 10 degrees tower occupancy for E particle > " + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "#eta";
    yleg = "#varphi";
    hHFminus10Occupancy[ihf] = MakeHisto(label,title,xleg,yleg,nHFtower_10,value_HFminus_eta_edge_10,nphi_10,value_HFminus_phi_edge_10);

    label = "hHFminus20Occupancy_" + TString::Format("%d",ihf+1);
    title = "HF minus 20 degrees tower occupancy for E particle > " + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "#eta";
    yleg = "#varphi";
    hHFminus20Occupancy[ihf] = MakeHisto(label,title,xleg,yleg,nHFtower_20,value_HFminus_eta_edge_20,nphi_20,value_HFminus_phi_edge_20);
  }

  cout<<endl<<endl;

  //-- HF plus correlation 
  double Nmoca_reco_sel_HFplus_correlation[nHFOccupancyCut];
  for(int i = 0; i < nHFOccupancyCut; i++) Nmoca_reco_sel_HFplus_correlation[i] = 0; 

  TH2D* hHFplus_correlation_energy[nHFOccupancyCut];
  TH2D* hHFplus_correlation_eta[nHFOccupancyCut];
  TH2D* hHFplus_correlation_phi[nHFOccupancyCut];

  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    label = "hHFplus_correlation_energy_" + TString::Format("%d",ihf+1);
    title = "HF plus energy correlation for E particle > "  + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "leading particle energy";
    yleg = "leading HF plus tower energy";
    hHFplus_correlation_energy[ihf] = MakeHisto(label,title,xleg,yleg,100,0,100,100,0,100);

    label = "hHFplus_correlation_eta_" + TString::Format("%d",ihf+1);									
    title = "HF plus eta correlation for E particle > "  + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";			
    xleg = "leading particle eta";
    yleg = "leading HF plus tower eta";
    hHFplus_correlation_eta[ihf] = MakeHisto(label,title,xleg,yleg,20,HFplus_eta_min,HFplus_eta_max,20,HFplus_eta_min,HFplus_eta_max);
    
    label = "hHFplus_correlation_phi_" + TString::Format("%d",ihf+1);							     
    title = "HF plus phi correlation for E particle > "  + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "leading particle phi";
    yleg = "leading HF plus tower phi";
    hHFplus_correlation_phi[ihf] = MakeHisto(label,title,xleg,yleg,50,-TMath::Pi(),TMath::Pi(),50,-TMath::Pi(),TMath::Pi());
  }

  cout<<endl<<endl;

  //-- HF minus correlation 
  double Nmoca_reco_sel_HFminus_correlation[nHFOccupancyCut];
  for(int i = 0; i < nHFOccupancyCut; i++) Nmoca_reco_sel_HFminus_correlation[i] = 0; 

  TH2D* hHFminus_correlation_energy[nHFOccupancyCut];
  TH2D* hHFminus_correlation_eta[nHFOccupancyCut];
  TH2D* hHFminus_correlation_phi[nHFOccupancyCut];

  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    label = "hHFminus_correlation_energy_" + TString::Format("%d",ihf+1);
    title = "HF minus energy correlation for E particle > "  + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "leading particle energy";
    yleg = "leading HF minus tower energy";
    hHFminus_correlation_energy[ihf] = MakeHisto(label,title,xleg,yleg,100,0,100,100,0,100);

    label = "hHFminus_correlation_eta_" + TString::Format("%d",ihf+1);									
    title = "HF minus eta correlation for E particle > "  + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";			
    xleg = "leading particle eta";
    yleg = "leading HF minus tower eta";
    hHFminus_correlation_eta[ihf] = MakeHisto(label,title,xleg,yleg,20,HFminus_eta_min,HFminus_eta_max,20,HFminus_eta_min,HFminus_eta_max);
    
    label = "hHFminus_correlation_phi_" + TString::Format("%d",ihf+1);							     
    title = "HF minus phi correlation for E particle > "  + TString::Format("%3.1f",EHFOccupancyCut[ihf]) + " GeV";
    xleg = "leading particle phi";
    yleg = "leading HF minus tower phi";
    hHFminus_correlation_phi[ihf] = MakeHisto(label,title,xleg,yleg,50,-TMath::Pi(),TMath::Pi(),50,-TMath::Pi(),TMath::Pi());
  }

  cout<<endl<<endl;

  //-- HF response
  TProfile *hHFplusResponseEgen = MakeProfile("hHFplusResponseEgen","HF plus response versus Egen","Egen","HF plus response",20,0,100);
  TProfile *hHFminusResponseEgen = MakeProfile("hHFminusResponseEgen","HF minus response versus Egen","Egen","HF minus response",20,0,100);
  
  TProfile *hHFplusResponseBiasEgen = MakeProfile("hHFplusResponseBiasEgen","HF plus response versus Egen","Egen","HF plus response",20,0,100);
  TProfile *hHFminusResponseBiasEgen = MakeProfile("hHFminusResponseBiasEgen","HF minus response versus Egen","Egen","HF minus response",20,0,100);
  
  cout<<endl<<endl;

  const int nbin_response = 10;
  TH1D *hHFResponse[nbin_response];
  TH1D *hHFplusResponse[nbin_response];
  TH1D *hHFminusResponse[nbin_response];

  for(int igen = 0; igen < nbin_response; ++igen) {
    double Emin = igen*100.0/nbin_response;
    double Emax = (igen+1)*100.0/nbin_response;
    
    label = "hHFResponse_" + TString::Format("%d",igen+1);
    title = "HF response for " + TString::Format("%3.1f",Emin) + " GeV" + " < E gen < "  + TString::Format("%3.1f",Emax) + " GeV";
    hHFResponse[igen] = MakeHisto(label,title,"response","N events",100,0,5); 

    label = "hHFplusResponse_" + TString::Format("%d",igen+1);
    title = "HF plus response for " + TString::Format("%3.1f",Emin) + " GeV" + " < E gen < "  + TString::Format("%3.1f",Emax) + " GeV";
    hHFplusResponse[igen] = MakeHisto(label,title,"response","N events",100,0,5);                                                       

    label = "hHFminusResponse_" + TString::Format("%d",igen+1);										
    title = "HF minus response for " + TString::Format("%3.1f",Emin) + " GeV" + " < E gen < "  + TString::Format("%3.1f",Emax) + " GeV";	  
    hHFminusResponse[igen] = MakeHisto(label,title,"response","N events",100,0,5);                                                       
  }

  cout<<endl<<endl;

  TH1D* hHFMeanResponse = MakeHisto("hHFMeanResponse","HF response versus Egen","Egen","HF  response",nbin_response,0,100);
  TH1D* hHFplusMeanResponse = MakeHisto("hHFplusMeanResponse","HF plus response versus Egen","Egen","HF plus response",nbin_response,0,100);
  TH1D* hHFminusMeanResponse = MakeHisto("hHFminusMeanResponse","HF minus response versus Egen","Egen","HF minus response",nbin_response,0,100);

  cout<<endl<<endl;
  
  //-- HF particle - tower matching													    
                                                                                                                                                  
  TH1D* hHFDeltaE = MakeHisto("hHFDeltaE","HF Delta E/E leading tower - matched tower","delta energy / energy","N events",100,0,100);  
  TH1D* hHFDeltaEta = MakeHisto("hHFDeltaEta","HF Delta Eta leading tower - matched tower","delta eta","N events",102,-5.05,5.15);	    
  TH1D* hHFDeltaPhi = MakeHisto("hHFDeltaPhi","HF Delta Phi leading tower - matched tower","delta phi","N events",100,-MY_PI,MY_PI);   

  //-- HF plus particle - tower matching

  TH1D* hHFplusDeltaE = MakeHisto("hHFplusDeltaE","HF plus Delta E/E leading tower - matched tower","delta energy / energy","N events",100,0,100);
  TH1D* hHFplusDeltaEta = MakeHisto("hHFplusDeltaEta","HF plus Delta Eta leading tower - matched tower","delta eta","N events",102,-5.05,5.15);
  TH1D* hHFplusDeltaPhi = MakeHisto("hHFplusDeltaPhi","HF plus Delta Phi leading tower - matched tower","delta phi","N events",100,-MY_PI,MY_PI);   

  //-- HF minus particle - tower matching

  TH1D* hHFminusDeltaE = MakeHisto("hHFminusDeltaE","HF minus Delta E/E leading tower - matched tower","delta energy / energy","N events",100,0,100);
  TH1D* hHFminusDeltaEta = MakeHisto("hHFminusDeltaEta","HF minus Delta Eta leading tower - matched tower","delta eta","N events",102,-5.05,5.15);
  TH1D* hHFminusDeltaPhi = MakeHisto("hHFminusDeltaPhi","HF minus Delta Phi leading tower - matched tower","delta phi","N events",100,-MY_PI,MY_PI);

  cout<<endl<<endl;

  //-- HF plus rechit

  TH1D *hHFplusRechitEnergyNoCut;											   
  label = "hHFplusRechitEnergyNoCut";											   
  title = "HF plus rechit energy no cut";										   
  hHFplusRechitEnergyNoCut = MakeHisto(label,title,"HF plus rechit energy","N HF plus rechits",1500,-4,26);		   
  
  TH1D *hnHFplusRechitNoCut;		  
  label = "hnHFplusRechitNoCut";
  title = "HF plus rechit multiplicity no cut";
  hnHFplusRechitNoCut = MakeHisto(label,title,"N HF plus rechits","N evts",2000,0,2000);                         

  cout<<endl;

  TH1D *hnHFplusRechit[nHFcut];		  
  TH1D *hHFplusRechitEtot[nHFcut];		
  TH1D *hHFplusRechitLeadingEnergy[nHFcut];	
                                     
  TH1D *hHFplusRechitEnergy[nHFcut];		
  TH1D *hHFplusRechitEt[nHFcut];
  TH1D *hHFplusRechitTime[nHFcut];  
  TH1D *hHFplusRechitEta[nHFcut];		
  TH1D *hHFplusRechitPhi[nHFcut];           
  TH2D *hHFplusRechitEtaPhi[nHFcut];
  
  TH1D *hHFplusRechitDepth[nHFcut];

  TH1D *hHFplusRechitiEtaDepth1[nHFcut];		 
  TH1D *hHFplusRechitiEtaDepth2[nHFcut];		 

  TH1D *hHFplusRechitiPhiDepth1[nHFcut];
  TH1D *hHFplusRechitiPhiDepth2[nHFcut];
  
  TH2D *hHFplusRechitiEtaiPhiDepth1[nHFcut];
  TH2D *hHFplusRechitiEtaiPhiDepth2[nHFcut];

  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    label = "hnHFplusRechit_" + TString::Format("%d",ihf+1);
    title = "HF plus rechit multiplicity for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hnHFplusRechit[ihf] = MakeHisto(label,title,"N HF plus rechits","N evts",2000,0,2000);                         

    label = "hHFplusRechitEtot_" + TString::Format("%d",ihf+1);						   
    title = "HF plus rechit total energy for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFplusRechitEtot[ihf] = MakeHisto(label,title,"HF plus rechit Etot","N evts",500,0,5000);
    
    label = "hHFplusRechitLeadingEnergy_" + TString::Format("%d",ihf+1);				     
    title = "HF plus rechit leading energy for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFplusRechitLeadingEnergy[ihf] = MakeHisto(label,title,"HF plus rechit leading energy","N evts",500,0,500);

    label = "hHFplusRechitEnergy_" + TString::Format("%d",ihf+1);						   
    title = "HF plus rechit energy for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFplusRechitEnergy[ihf] = MakeHisto(label,title,"HF plus rechit energy","N HF plus rechits",500,0,500);

    label = "hHFplusRechitEt_" + TString::Format("%d",ihf+1);						   
    title = "HF plus rechit Et for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFplusRechitEt[ihf] = MakeHisto(label,title,"HF plus rechit Et","N HF plus rechits",200,0,20);

    label = "hHFplusRechitTime_" + TString::Format("%d",ihf+1);						    
    title = "HF plus rechit time for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	    
    hHFplusRechitTime[ihf] = MakeHisto(label,title,"HF plus rechit time","N HF plus rechits",280,-11,17);

    label = "hHFplusRechitEta_" + TString::Format("%d",ihf+1);				     
    title = "HF plus rechit eta for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFplusRechitEta[ihf] = MakeHisto(label,title,"HF plus rechit eta","N HF plus rechits",120,-6,6); 

    label = "hHFplusRechitPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF plus rechit phi for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFplusRechitPhi[ihf] = MakeHisto(label,title,"HF plus rechit phi","N HF plus rechits",100,-TMath::Pi(),TMath::Pi());

    label = "hHFplusRechitEtaPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF plus rechit eta-phi for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFplusRechitEtaPhi[ihf] = MakeHisto(label,title,"HF plus rechit eta","HF plus rechit phi",120,-6,6,100,-TMath::Pi(),TMath::Pi());

    label = "hHFplusRechitDepth_" + TString::Format("%d",ihf+1);					     
    title = "HF plus rechit depth for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	     
    hHFplusRechitDepth[ihf] = MakeHisto(label,title,"HF plus rechit depth","N HF plus rechits",2,0.5,2.5);   

    label = "hHFplusRechitiEtaDepth1_" + TString::Format("%d",ihf+1);
    title = "HF plus rechit ieta depth 1 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFplusRechitiEtaDepth1[ihf] = MakeHisto(label,title,"HF plus rechit ieta depth 1","N HF plus rechits",15,27.5,42.5);

    label = "hHFplusRechitiEtaDepth2_" + TString::Format("%d",ihf+1);							 
    title = "HF plus rechit ieta depth 2 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 
    hHFplusRechitiEtaDepth2[ihf] = MakeHisto(label,title,"HF plus rechit ieta depth 2","N HF plus rechits",15,27.5,42.5);   

    label = "hHFplusRechitiPhiDepth1_" + TString::Format("%d",ihf+1);
    title = "HF plus rechit iphi depth 1 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFplusRechitiPhiDepth1[ihf] = MakeHisto(label,title,"HF plus rechit iphi depth 1","N HF plus rechits",71,1,72);

    label = "hHFplusRechitiPhiDepth2_" + TString::Format("%d",ihf+1);							 
    title = "HF plus rechit iphi depth 2 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 
    hHFplusRechitiPhiDepth2[ihf] = MakeHisto(label,title,"HF plus rechit iphi depth 2","N HF plus rechits",71,1,72);  

    label = "hHFplusRechitiEtaiPhiDepth1_" + TString::Format("%d",ihf+1);							 
    title = "HF plus rechit ieta-iphi depth 1 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 
    hHFplusRechitiEtaiPhiDepth1[ihf] = MakeHisto(label,title,"HF plus rechit ieta depth 1","HF plus rechit iphi depth 1",15,27.5,42.5,71,1,72);
   
    label = "hHFplusRechitiEtaiPhiDepth2_" + TString::Format("%d",ihf+1);							 		       
    title = "HF plus rechit ieta-iphi depth 2 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 			  
    hHFplusRechitiEtaiPhiDepth2[ihf] = MakeHisto(label,title,"HF plus rechit ieta depth 2","HF plus rechit iphi depth 2",15,27.5,42.5,71,1,72);   
    
    cout<<endl<<endl;
  }
  
  //-- HF minus rechit

  TH1D *hHFminusRechitEnergyNoCut;											   
  label = "hHFminusRechitEnergyNoCut";											   
  title = "HF minus rechit energy no cut";										   
  hHFminusRechitEnergyNoCut = MakeHisto(label,title,"HF minus rechit energy","N HF minus rechits",1500,-4,26);		   

  TH1D *hnHFminusRechitNoCut;		  									 
  label = "hnHFminusRechitNoCut";										 
  title = "HF minus rechit multiplicity no cut";									 
  hnHFminusRechitNoCut = MakeHisto(label,title,"N HF minus rechits","N evts",2000,0,2000);                         

  cout<<endl;

  TH1D *hnHFminusRechit[nHFcut];		  
  TH1D *hHFminusRechitEtot[nHFcut];		
  TH1D *hHFminusRechitLeadingEnergy[nHFcut];	
                                     
  TH1D *hHFminusRechitEnergy[nHFcut];		
  TH1D *hHFminusRechitEt[nHFcut];
  TH1D *hHFminusRechitTime[nHFcut];  
  TH1D *hHFminusRechitEta[nHFcut];		
  TH1D *hHFminusRechitPhi[nHFcut];           
  TH2D *hHFminusRechitEtaPhi[nHFcut];
  
  TH1D *hHFminusRechitDepth[nHFcut];

  TH1D *hHFminusRechitiEtaDepth1[nHFcut];		 
  TH1D *hHFminusRechitiEtaDepth2[nHFcut];		 

  TH1D *hHFminusRechitiPhiDepth1[nHFcut];
  TH1D *hHFminusRechitiPhiDepth2[nHFcut];
  
  TH2D *hHFminusRechitiEtaiPhiDepth1[nHFcut];
  TH2D *hHFminusRechitiEtaiPhiDepth2[nHFcut];

  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    label = "hnHFminusRechit_" + TString::Format("%d",ihf+1);
    title = "HF minus rechit multiplicity for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hnHFminusRechit[ihf] = MakeHisto(label,title,"N HF minus rechits","N evts",2000,0,2000);                         

    label = "hHFminusRechitEtot_" + TString::Format("%d",ihf+1);						   
    title = "HF minus rechit total energy for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFminusRechitEtot[ihf] = MakeHisto(label,title,"HF minus rechit Etot","N evts",500,0,5000);
    
    label = "hHFminusRechitLeadingEnergy_" + TString::Format("%d",ihf+1);				     
    title = "HF minus rechit leading energy for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFminusRechitLeadingEnergy[ihf] = MakeHisto(label,title,"HF minus rechit leading energy","N evts",500,0,500);

    label = "hHFminusRechitEnergy_" + TString::Format("%d",ihf+1);						   
    title = "HF minus rechit energy for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFminusRechitEnergy[ihf] = MakeHisto(label,title,"HF minus rechit energy","N HF minus rechits",500,0,500);

    label = "hHFminusRechitEt_" + TString::Format("%d",ihf+1);						   
    title = "HF minus rechit Et for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFminusRechitEt[ihf] = MakeHisto(label,title,"HF minus rechit Et","N HF minus rechits",200,0,20);

    label = "hHFminusRechitTime_" + TString::Format("%d",ihf+1);						    
    title = "HF minus rechit time for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	    
    hHFminusRechitTime[ihf] = MakeHisto(label,title,"HF minus rechit time","N HF minus rechits",280,-11,17);

    label = "hHFminusRechitEta_" + TString::Format("%d",ihf+1);				     
    title = "HF minus rechit eta for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFminusRechitEta[ihf] = MakeHisto(label,title,"HF minus rechit eta","N HF minus rechits",120,-6,6); 

    label = "hHFminusRechitPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF minus rechit phi for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFminusRechitPhi[ihf] = MakeHisto(label,title,"HF minus rechit phi","N HF minus rechits",100,-TMath::Pi(),TMath::Pi());

    label = "hHFminusRechitEtaPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF minus rechit eta-phi for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFminusRechitEtaPhi[ihf] = MakeHisto(label,title,"HF minus rechit eta","HF minus rechit phi",120,-6,6,100,-TMath::Pi(),TMath::Pi());

    label = "hHFminusRechitDepth_" + TString::Format("%d",ihf+1);					     
    title = "HF minus rechit depth for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	     
    hHFminusRechitDepth[ihf] = MakeHisto(label,title,"HF minus rechit depth","N HF minus rechits",2,0.5,2.5);   

    label = "hHFminusRechitiEtaDepth1_" + TString::Format("%d",ihf+1);
    title = "HF minus rechit ieta depth 1 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFminusRechitiEtaDepth1[ihf] = MakeHisto(label,title,"HF minus rechit ieta depth 1","N HF minus rechits",15,-42.5,-27.5);

    label = "hHFminusRechitiEtaDepth2_" + TString::Format("%d",ihf+1);							 
    title = "HF minus rechit ieta depth 2 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 
    hHFminusRechitiEtaDepth2[ihf] = MakeHisto(label,title,"HF minus rechit ieta depth 2","N HF minus rechits",15,-42.5,-27.5); 

    label = "hHFminusRechitiPhiDepth1_" + TString::Format("%d",ihf+1);
    title = "HF minus rechit iphi depth 1 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFminusRechitiPhiDepth1[ihf] = MakeHisto(label,title,"HF minus rechit iphi depth 1","N HF minus rechits",71,1,72);

    label = "hHFminusRechitiPhiDepth2_" + TString::Format("%d",ihf+1);							 
    title = "HF minus rechit iphi depth 2 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 
    hHFminusRechitiPhiDepth2[ihf] = MakeHisto(label,title,"HF minus rechit iphi depth 2","N HF minus rechits",71,1,72);   

    label = "hHFminusRechitiEtaiPhiDepth1_" + TString::Format("%d",ihf+1);							 
    title = "HF minus rechit ieta-iphi depth 1 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 
    hHFminusRechitiEtaiPhiDepth1[ihf] = MakeHisto(label,title,"HF minus rechit ieta depth 1","HF minus rechit iphi depth 1",
						  15,-42.5,-27.5,71,1,72);   

    label = "hHFminusRechitiEtaiPhiDepth2_" + TString::Format("%d",ihf+1);							 		       
    title = "HF minus rechit ieta-iphi depth 2 for E rechit > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		 			
    hHFminusRechitiEtaiPhiDepth2[ihf] = MakeHisto(label,title,"HF minus rechit ieta depth 2","HF minus rechit iphi depth 2",
						  15,-42.5,-27.5,71,1,72);   
    
    cout<<endl<<endl;
  }
  
  //-- HF tower

  TH1D *hnHFTower[nHFcut];
  TH1D *hHFTowerEtot[nHFcut];
  TH1D *hHFTowerLeadingEnergy[nHFcut];

  TH1D *hHFTowerEnergy[nHFcut];
  TH1D *hHFTowerEta[nHFcut];
  TH1D *hHFTowerPhi[nHFcut];            

  TH2D *hHFTowerEtaPhi[nHFcut];

  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    label = "hnHFTower_" + TString::Format("%d",ihf+1);
    title = "HF tower multiplicity for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hnHFTower[ihf] = MakeHisto(label,title,"N HF towers","N evts",500,0,500);                         

    label = "hHFTowerEtot_" + TString::Format("%d",ihf+1);						   
    title = "HF tower total energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFTowerEtot[ihf] = MakeHisto(label,title,"HF tower Etot","N evts",500,0,5000);
    
    label = "hHFTowerLeadingEnergy_" + TString::Format("%d",ihf+1);				     
    title = "HF tower leading energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFTowerLeadingEnergy[ihf] = MakeHisto(label,title,"HF tower leading energy","N evts",500,0,500);

    label = "hHFTowerEnergy_" + TString::Format("%d",ihf+1);						   
    title = "HF tower energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFTowerEnergy[ihf] = MakeHisto(label,title,"HF tower energy","N HF towers",500,0,500);

    label = "hHFTowerEta_" + TString::Format("%d",ihf+1);				     
    title = "HF tower eta for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFTowerEta[ihf] = MakeHisto(label,title,"HF tower eta","N HF towers",120,-6,6); 

    label = "hHFTowerPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF tower phi for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFTowerPhi[ihf] = MakeHisto(label,title,"HF tower phi","N HF towers",100,-TMath::Pi(),TMath::Pi());

    label = "hHFTowerEtaPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF tower eta-phi for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFTowerEtaPhi[ihf] = MakeHisto(label,title,"HF tower eta","HF tower phi",120,-6,6,100,-TMath::Pi(),TMath::Pi());

    cout<<endl<<endl;
  }
                                                                                                                                       
  //-- HF plus tower												   
  
  TH1D *hHFplusTowerEnergyNoCut;
  label = "hHFplusTowerEnergyNoCut";
  title = "HF plus tower energy no cut";
  hHFplusTowerEnergyNoCut = MakeHisto(label,title,"HF plus tower energy","N HF plus towers",1500,-4,26);		   

  TH1D *hnHFplusTowerNoCut;										   
  label = "hnHFplusTowerNoCut";
  title = "HF plus tower multiplicity no cut";
  hnHFplusTowerNoCut = MakeHisto(label,title,"N HF plus towers","N evts",500,0,500);				   

  cout<<endl;

  TH1D *hnHFplusTower[nHFcut];										   
  TH1D *hHFplusTowerEtot[nHFcut];										   
  TH1D *hHFplusTowerLeadingEnergy[nHFcut];									   
  
  TH1D *hHFplusTowerEnergy[nHFcut];										   
  TH1D *hHFplusTowerEta[nHFcut];										   
  TH1D *hHFplusTowerPhi[nHFcut];

  TH2D *hHFplusTowerEtaPhi[nHFcut];
  
  for(int ihf = 0; ihf < nHFcut; ++ihf) {								   
    label = "hnHFplusTower_" + TString::Format("%d",ihf+1);							   
    title = "HF plus tower multiplicity for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	   
    hnHFplusTower[ihf] = MakeHisto(label,title,"N HF plus towers","N evts",500,0,500);				   
    
    label = "hHFplusTowerEtot_" + TString::Format("%d",ihf+1);						   
    title = "HF plus tower total energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	   
    hHFplusTowerEtot[ihf] = MakeHisto(label,title,"HF plus tower Etot","N evts",500,0,5000);			   
    
    label = "hHFplusTowerLeadingEnergy_" + TString::Format("%d",ihf+1);				     	   
    title = "HF plus tower leading energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	   
    hHFplusTowerLeadingEnergy[ihf] = MakeHisto(label,title,"HF plus tower leading energy","N evts",500,0,500);	   
    
    label = "hHFplusTowerEnergy_" + TString::Format("%d",ihf+1);						   
    title = "HF plus tower energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		   
    hHFplusTowerEnergy[ihf] = MakeHisto(label,title,"HF plus tower energy","N HF plus towers",500,0,500);		   
    
    label = "hHFplusTowerEta_" + TString::Format("%d",ihf+1);				     		   
    title = "HF plus tower eta for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    		   
    hHFplusTowerEta[ihf] = MakeHisto(label,title,"HF plus tower eta","N HF plus towers",120,-6,6); 			   
                                                                                                         
    label = "hHFplusTowerPhi_" + TString::Format("%d",ihf+1);				          	   
    title = "HF plus tower phi for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    		   
    hHFplusTowerPhi[ihf] = MakeHisto(label,title,"HF plus tower phi","N HF plus towers",100,-TMath::Pi(),TMath::Pi());   

    label = "hHFplusTowerEtaPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF plus tower eta-phi for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFplusTowerEtaPhi[ihf] = MakeHisto(label,title,"HF plus tower eta","HF plus tower phi",120,-6,6,100,-TMath::Pi(),TMath::Pi());

    cout<<endl<<endl;
  }                                                                                                                     

  //-- HF minus tower	

  TH1D *hHFminusTowerEnergyNoCut;											   
  label = "hHFminusTowerEnergyNoCut";											   
  title = "HF minus tower energy no cut";										   
  hHFminusTowerEnergyNoCut = MakeHisto(label,title,"HF minus tower energy","N HF minus towers",1500,-4,26);		 

  TH1D *hnHFminusTowerNoCut;										   	   
  label = "hnHFminusTowerNoCut";											   
  title = "HF minus tower multiplicity no cut";									   
  hnHFminusTowerNoCut = MakeHisto(label,title,"N HF minus towers","N evts",500,0,500);				   

  cout<<endl;

  TH1D *hnHFminusTower[nHFcut];										   		
  TH1D *hHFminusTowerEtot[nHFcut];										   	
  TH1D *hHFminusTowerLeadingEnergy[nHFcut];									   	
															
  TH1D *hHFminusTowerEnergy[nHFcut];										   	
  TH1D *hHFminusTowerEta[nHFcut];										   	
  TH1D *hHFminusTowerPhi[nHFcut];	

  TH2D *hHFminusTowerEtaPhi[nHFcut];
															
  for(int ihf = 0; ihf < nHFcut; ++ihf) {								   		
    label = "hnHFminusTower_" + TString::Format("%d",ihf+1);							   	
    title = "HF minus tower multiplicity for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	   	
    hnHFminusTower[ihf] = MakeHisto(label,title,"N HF minus towers","N evts",500,0,500);				   	
    
    label = "hHFminusTowerEtot_" + TString::Format("%d",ihf+1);						   		
    title = "HF minus tower total energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	   	
    hHFminusTowerEtot[ihf] = MakeHisto(label,title,"HF minus tower Etot","N evts",500,0,5000);			   	
    
    label = "hHFminusTowerLeadingEnergy_" + TString::Format("%d",ihf+1);				     	   		
    title = "HF minus tower leading energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";	   	
    hHFminusTowerLeadingEnergy[ihf] = MakeHisto(label,title,"HF minus tower leading energy","N evts",500,0,500);	   	
    
    label = "hHFminusTowerEnergy_" + TString::Format("%d",ihf+1);						   	
    title = "HF minus tower energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";		   	
    hHFminusTowerEnergy[ihf] = MakeHisto(label,title,"HF minus tower energy","N HF minus towers",500,0,500);		
    
    label = "hHFminusTowerEta_" + TString::Format("%d",ihf+1);				     		   		
    title = "HF minus tower eta for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    		   	
    hHFminusTowerEta[ihf] = MakeHisto(label,title,"HF minus tower eta","N HF minus towers",120,-6,6); 			
    
    label = "hHFminusTowerPhi_" + TString::Format("%d",ihf+1);				          	   		
    title = "HF minus tower phi for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    		   	
    hHFminusTowerPhi[ihf] = MakeHisto(label,title,"HF minus tower phi","N HF minus towers",100,-TMath::Pi(),TMath::Pi());  

    label = "hHFminusTowerEtaPhi_" + TString::Format("%d",ihf+1);				          
    title = "HF minus tower eta-phi for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";    
    hHFminusTowerEtaPhi[ihf] = MakeHisto(label,title,"HF minus tower eta","HF minus tower phi",120,-6,6,100,-TMath::Pi(),TMath::Pi());

    cout<<endl<<endl;
  }                                                                                                                     

  //-- correlation coefficient
  
  const double delta_eta_0 = 0.5;
  const int n_delta_eta = 9;
  
  double eta_plus_min[n_delta_eta];
  double eta_plus_max[n_delta_eta];
  
  double eta_minus_min[n_delta_eta];
  double eta_minus_max[n_delta_eta];

  TH1D *hnFwd[n_delta_eta];
  TH1D *hnBwd[n_delta_eta];
  TH1D *hnFwdBwd[n_delta_eta];
  
  for(int i = 0; i < n_delta_eta; ++i) {
    double delta_eta = i*delta_eta_0;
    
    eta_plus_min[i] = 0.5*delta_eta;
    eta_plus_max[i] = 0.5*delta_eta + delta_eta_0;
    //-- cout<<"forward interval: ["<<eta_plus_min[i]<<","<<eta_plus_max[i]<<"]"<<endl; 
    
    eta_minus_min[i] = -0.5*delta_eta - delta_eta_0;
    eta_minus_max[i] = -0.5*delta_eta;
    //-- cout<<"backward interval: ["<<eta_minus_min[i]<<","<<eta_minus_max[i]<<"]"<<endl;

    label = "hnFwd_" + TString::Format("%d",i);
    title = "n fwd for " + TString::Format("%3.2f",eta_plus_min[i]) + " < eta < " + TString::Format("%3.2f",eta_plus_max[i]);
    hnFwd[i] = MakeHisto(label,title,"N fwd tracks","N evts",100,0,100);                                                              
    
    label = "hnBwd_" + TString::Format("%d",i);											      
    title = "n bwd for " + TString::Format("%3.2f",eta_minus_min[i]) + " < eta < " + TString::Format("%3.2f",eta_minus_max[i]);
    hnBwd[i] = MakeHisto(label,title,"N bwd tracks","N evts",100,0,100);                                                                 

    label = "hnFwdBwd_" + TString::Format("%d",i);											     
    title = "n fwd * n bwd for " + TString::Format("%3.2f",eta_plus_min[i]) + " < eta < " + TString::Format("%3.2f",eta_plus_max[i]); 
    hnFwdBwd[i] = MakeHisto(label,title,"N fwd tracks * N bwd tracks","N evts",100,0,2000);                                                     
            
    //-- getchar();
  }

  cout<<endl;

  //-- matched charged particle distributions - exactly one good vertex
  TH1D *hMChPt = MakeHisto("hMChPt","charged particle pT","charged particle pT","N charged particles",150,0,15);			    
  TH1D *hMChEta = MakeHisto("hMChEta","charged particle eta","charged particle #eta","N charged particles",120,-3,3);		      
  TH1D *hMChPhi = MakeHisto("hMChPhi","charged particle phi","charged particle #phi","N charged particles",100,-TMath::Pi(),TMath::Pi());
  
  TH2D *hMChMTPt = MakeHisto("hMChMTPt","track pT versus charged particle pT","charged particle pT","track pT",150,0,15,150,0,15);
  TH2D *hMChMTEta = MakeHisto("hMChMTEta","track eta versus charged particle eta","charged particle eta","track eta",120,-3,3,120,-3,3);
  TH2D *hMChMTPhi = MakeHisto("hMChMTPhi","track phi versus charged particle phi","charged particle phi","track phi",100,-MY_PI,MY_PI,100,-MY_PI,MY_PI);

  cout<<endl;

  //-- track distributions exactly one good vertex
  TH1D *hTrackNum = MakeHisto("hTrackNum","track multiplicity","N tracks","N evts",150,0,150);
  TH1D *hTrackPt = MakeHisto("hTrackPt","track pT","track pT","N tracks",150,0,15);
  TH1D *hTrackEta = MakeHisto("hTrackEta","track eta","track #eta","N tracks",120,-3,3);

  TH1D *hTrackPhi = MakeHisto("hTrackPhi","track phi","track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());
  TH1D *hTrackPhiCentral = MakeHisto("hTrackPhiCentral","track phi |eta| < 1","track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());
  TH1D *hTrackPhiFwd = MakeHisto("hTrackPhiFwd","track phi |eta| > 1","track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());

  TH2D* hTrackEtaPhi = MakeHisto("hTrackEtaPhi","track eta - phi","track #eta","track #phi",240,-3,3,200,-TMath::Pi(),TMath::Pi());
  TH2D* hTrackEtaPhiPixel = MakeHisto("hTrackEtaPhiPixel","track eta - phi pixel","track #eta","track #phi",240,-3,3,200,-TMath::Pi(),TMath::Pi());
  TH2D* hTrackEtaPhiStrip = MakeHisto("hTrackEtaPhiStrip","track eta - phi strip","track #eta","track #phi",240,-3,3,200,-TMath::Pi(),TMath::Pi());

  TH1D *hTrackQuality = MakeHisto("hTrackQuality","track quality","track quality","N tracks",3,-0.5,2.5);

  const int neta_edge = 5;
  double eta_edge[neta_edge] = {0,0.6,1.2,1.8,2.4};
  
  TH1D *hTrackPtEta[neta_edge-1];
  
  for(int ieta = 0; ieta < neta_edge-1; ieta++) {
    label = "hTrackPtEta_" + TString::Format("%d",ieta+1);
    title = " track pT for " + TString::Format("%3.1f",eta_edge[ieta]) + " < |eta| < " + TString::Format("%3.1f",eta_edge[ieta+1]);
    hTrackPtEta[ieta] = MakeHisto(label,title,"track pT","N tracks",150,0,15);
  }

  TH1D *hTrackPtError = MakeHisto("hTrackPtError","track pT error/pT","track pT error/pT","N tracks",100,0,10);
  TH1D *hTrackEtaError = MakeHisto("hTrackEtaError","track eta error","track #eta error", "N tracks", 100,0,0.02);
  TH1D *hTrackPhiError = MakeHisto("hTrackPhiError","track phi error","track #phi error","N tracks",100,0,0.02);
  
  TProfile *hTrackPtErrorPt = MakeProfile("hTrackPtErrorPt","track pT error/pT versus pT","track pT","track pT error/pT",100,0,5);
  TProfile *hTrackPtErrorEta = MakeProfile("hTrackPtErrorEta","track pT error/pT versus #eta","track #eta","track pT error/pT",100,-2.5,2.5);
  TProfile *hTrackPtErrorNhit = MakeProfile("hTrackPtErrorNhit","track pT error/pT versus Nhit","track Nhit","track pT error/pT",51,-0.5,50.5);

  TH1D *hTrackdxy = MakeHisto("hTrackdxy","track dxy","track dxy/sigmaxy","N tracks",100,0,10);
  TH1D *hTrackdz = MakeHisto("hTrackdz","track dz","track dz/sigmaz","N tracks",100,0,10);

  TProfile *hTrackdxyEta = MakeProfile("hTrackdxyEta","track dxy/sigmaxy versus eta","track eta","track dxy/sigmaxy",120,-3,3);
  TProfile *hTrackdzEta = MakeProfile("hTrackdzEta","track dz/sigmaz versus eta","track eta","track dz/sigmaz",120,-3,3);       

  TProfile *hTrackdxyNhit = MakeProfile("hTrackdxyNhit","track dxy/sigmaxy versus Nhit","track Nhit","track dxy/sigmaxy",51,-0.5,50.5);	
  TProfile *hTrackdzNhit = MakeProfile("hTrackdzNhit","track dz/sigmaz versus Nhit","track Nhit","track dz/sigmaz",51,-0.5,50.5);       
  
  TH2D *hTrackdxyEta2D = MakeHisto("hTrackdxyEta2D","track dxy/sigmaxy versus eta","track eta","track dxy/sigmaxy",240,-3,3,200,-10,10);
  TH2D *hTrackdzEta2D = MakeHisto("hTrackdzEta2D","track dz/sigmaz versus eta","track eta","track dz/sigmaz",240,-3,3,200,-10,10);

  TH2D *hTrackdxyNum2D = MakeHisto("hTrackdxyNum2D","track dxy/sigmaxy versus N tracks","N tracks","track dxy/sigmaxy",200,0,100,200,-10,10);
  TH2D *hTrackdzNum2D = MakeHisto("hTrackdzNum2D","track dz/sigmaz versus N tracks","N tracks","track dz/sigmaz",200,0,100,200,-10,10);

  TH1D *hTrackchi2 = MakeHisto("hTrackchi2","track chi2","track chi2","N tracks",100,0,100);
  TH1D *hTrackndof = MakeHisto("hTrackndof","track ndof","track ndof","N tracks",100,0,100);
  TH1D *hTrackchi2ndof = MakeHisto("hTrackchi2ndof","track chi2ndof","track chi2ndof","N tracks",100,0,10);

  TH1D *hTrackValidHit = MakeHisto("hTrackValidHit","track ValidHit","track ValidHit","N tracks",101,-0.5,100.5);
  TH1D *hTrackValidPixelHit = MakeHisto("hTrackValidPixelHit","track ValidPixelHit","track ValidPixelHit","N tracks",101,-0.5,100.5);
  TH1D *hTrackValidStripHit = MakeHisto("hTrackValidStripHit","track ValidStripHit","track ValidStripHit","N tracks",101,-0.5,100.5);
 
  TH1D *hTrackWeight = MakeHisto("hTrackWeight","track weight","track weight","N tracks",100,0,1);

  TProfile *hTrackValidHitEta = MakeProfile("hTrackValidHitEta","track ValidHit versus eta","eta","track ValidHit",120,-3,3);
  TProfile *hTrackValidPixelHitEta = MakeProfile("hTrackValidPixelHitEta","track ValidPixelHit versus eta","eta","track ValidPixelHit",120,-3,3);
  TProfile *hTrackValidStripHitEta = MakeProfile("hTrackValidStripHitEta","track ValidStripHit versus eta","eta","track ValidStripHit",120,-3,3);

  TH1D *hLeadingTrackPt = MakeHisto("hLeadingTrackPt","leading track pT","leading track pT","N tracks",150,0,15);
  TH1D *hLeadingTrackEta = MakeHisto("hLeadingTrackEta","leading track eta","leading track #eta","N tracks",120,-3,3);
  TH1D *hLeadingTrackPhi = MakeHisto("hLeadingTrackPhi","leading track phi","leading track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());
  cout<<endl;

  //-- Track Hit 
  const int detsize = 6;

  TString DetName[detsize];
  DetName[0] = "PXB";
  DetName[1] = "PXF";
  DetName[2] = "TIB";
  DetName[3] = "TID";
  DetName[4] = "TOB";
  DetName[5] = "TEC";

  TString SubdetName[detsize];	
  SubdetName[0] = "layer";	
  SubdetName[1] = "disk";	
  SubdetName[2] = "layer";	
  SubdetName[3] = "wheel";	
  SubdetName[4] = "layer";	
  SubdetName[5] = "wheel";   

  int MaxLayer[detsize];
  for (int idet = 0; idet < detsize; idet++) MaxLayer[idet] = 0;   

  //-- all tracker
  TH2D* hHitEtaPhiAll = MakeHisto("hHitEtaPhiAll","Hit phi versus eta","eta","phi",240,-3,3,200,-TMath::Pi(),TMath::Pi());  
  TH2D* hHitZRAll = MakeHisto("hHitZRAll","Hit R versus z","z","R",1200,-300,300,1200,0,150);                               
  TH2D* hHitXYAll = MakeHisto("hHitXYAll","Hit y versus x","x","y",1200,-150,150,1200,-150,150);                        

  cout<<endl;

  //-- each detector
  TH2D *hHitEtaPhiDet[detsize];
  TH2D *hHitZRDet[detsize];
  TH2D *hHitXYDet[detsize];    
  
  double bRZ[detsize][4] = {{-30,30,3,12},{-60,60,5,15},{-80,80,20,55},{-150,150,25,50},{-120,120,50,120},{-300,300,20,110}};
  double nRZ[detsize][2] = {{600,90},{1200,100},{1600,350},{3000,250},{2400,700},{6000,900}};

  double bXY[detsize][4] = {{-15,15,-15,15},{-20,20,-20,20},{-60,60,-60,60},{-50,50,-50,50},{-125,125,-125,125},{-110,110,-110,110}};
  double nXY[detsize][2] = {{300,300},{400,400},{1200,1200},{1000,1000},{2500,2500},{2200,2200}};

  for (int idet = 0; idet < detsize; idet++){
    label = "hHitEtaPhi" + DetName[idet];
    title = "Hit phi versus eta - " + DetName[idet];
    hHitEtaPhiDet[idet] = MakeHisto(label,title,"eta","phi",240,-3,3,200,-TMath::Pi(),TMath::Pi()); 

    label = "hHitZR" + DetName[idet];
    title = "Hit R versus z - " + DetName[idet];
    hHitZRDet[idet] = MakeHisto(label,title,"z","R",nRZ[idet][0],bRZ[idet][0],bRZ[idet][1],nRZ[idet][1],bRZ[idet][2],bRZ[idet][3]);

    label = "hHitXY" + DetName[idet];
    title = "Hit y versus x - " + DetName[idet];
    hHitXYDet[idet] = MakeHisto(label,title,"x","y",nXY[idet][0],bXY[idet][0],bXY[idet][1],nXY[idet][1],bXY[idet][2],bXY[idet][3]);
  }          
  
  cout<<endl;

  //-- each layer
  int nlayer[detsize];
  nlayer[0] = 3;
  nlayer[1] = 2;
  nlayer[2] = 4;
  nlayer[3] = 3;
  nlayer[4] = 6;
  nlayer[5] = 9;

  int nlayermax = 9;
  TH2D *hHitEtaPhiLayer[detsize][nlayermax]; 	
  TH2D *hHitZRLayer[detsize][nlayermax]; 	
  TH2D *hHitXYLayer[detsize][nlayermax];      
  
  for (int idet = 0; idet < detsize; idet++){
    for (int ilayer = 0; ilayer < nlayer[idet]; ilayer++) {
      label = "hHitEtaPhi" + DetName[idet] + "_" + SubdetName[idet] + "_" + TString::Format("%d",ilayer+1);					 
      title = "Hit phi versus eta - " + DetName[idet] + " - " + SubdetName[idet] + " - " + TString::Format("%d",ilayer+1);    
      hHitEtaPhiLayer[idet][ilayer] = MakeHisto(label,title,"eta","phi",240,-3,3,200,-TMath::Pi(),TMath::Pi()); 	
      													
      label = "hHitZR" + DetName[idet] + "_" + SubdetName[idet] + "_" + TString::Format("%d",ilayer+1);
      title = "Hit R versus z - " + DetName[idet] + " - " + SubdetName[idet] + " - " + TString::Format("%d",ilayer+1); 				   
      hHitZRLayer[idet][ilayer] = MakeHisto(label,title,"z","R",1200,-300,300,1200,0,150);         			
      													
      label = "hHitXY" + DetName[idet] + "_" + SubdetName[idet] + "_" + TString::Format("%d",ilayer+1);
      title = "Hit y versus x - " + DetName[idet] + "_" + SubdetName[idet] + "_" + TString::Format("%d",ilayer+1);			  
      hHitXYLayer[idet][ilayer] = MakeHisto(label,title,"x","y",1200,-150,150,1200,-150,150);                   
    }
  }          
  
  cout<<endl;

  //-- vertex distributions exactly one good vertex
  TH1D *hVx = MakeHisto("hVx","vertex x","vertex x [cm]","N evts",200,-0.21,0.19);
  TH1D *hVxBS = MakeHisto("hVxBS","vertex x BS","vertex x BS [cm]","N evts",800,-0.2,0.2);
  TH1D *hVy = MakeHisto("hVy","vertex y","vertex y [cm]","N evts",200,-0.21,0.19);
  TH1D *hVyBS = MakeHisto("hVyBS","vertex y BS","vertex y BS [cm]","N evts",800,-0.2,0.2);
  TH2D *hVxy = MakeHisto("hVxy","vertex y versus x","vertex x [cm]","vertex y [cm]",200,-0.21,0.19,200,-0.21,0.19);
  TH1D *hVz = MakeHisto("hVz","vertex z","vertex z [cm]","N evts",210,-21.5,20.5);
  TH1D *hVrho = MakeHisto("hVrho","vertex rho","vertex #rho [cm]","N evts",2000,0,1);
  TH1D *hVrhoBS = MakeHisto("hVrhoBS","vertex rho BS","vertex #rho BS [cm]","N evts",10000,0,1);
  TH1D *hVchi2 = MakeHisto("hVchi2","vertex chi2","vertex #chi^{2}","N evts",500,0,100);
  TH1D *hVndof = MakeHisto("hVndof","vertex ndof","vertex ndof","N evts",100,0,100);
  TH1D *hVchi2ndof = MakeHisto("hVchi2ndof","vertex chi2/ndof","vertex #chi^{2}/ndof","N evts",200,0,10); 
  TH1D *hVzError = MakeHisto("hVzError","vertex z error","vertex z error [cm]","N evts",200,0,0.2);
  TH2D *hVzErrorVz = MakeHisto("hVzErrorVz","vertex z error versus vertex z","vertex z [cm]","vertex z error [cm]",100,-15,15,100,0,0.2);
  TH2D *hVzErrorTrackNum = MakeHisto("hVzErrorTrackNum","vertex z error versus track multiplicity","N tracks","vertex z error [cm]",100,0,100,100,0,0.2);
  cout<<endl;

  //-- beam spot distributions
  TH1D *hBSx = MakeHisto("hBSx","BS x","BS x [cm]","N evts",1000,-0.21,0.19);
  TH1D *hBSy = MakeHisto("hBSy","BS y","BS y [cm]","N evts",1000,-0.21,0.19);
  TH1D *hBSz = MakeHisto("hBSz","BS z","BS z [cm]","N evts",1000,-20.05,19.95);

  //-- all vertices
  TH1D *hVNum = MakeHisto("hVNum","vertex multiplicity","vertex multiplicity","N events",10,0.5,10.5);
  TH1D *hVNumGood = MakeHisto("hVNumGood","good vertex multiplicity","good vertex multiplicity","N events",10,0.5,10.5);
  TH2D *hVNum2D = MakeHisto("hVNum2D","n good vertex versus n vertex","n vertex","n good vertex",10,0.5,10.5,10,0.5,10.5);
  
  //-- primary vertices
  TH1D *hPVx = MakeHisto("hPVx","primary vertex x","primary vertex x [cm]","N evts",200,-0.21,0.19);
  TH1D *hPVy = MakeHisto("hPVy","primary vertex y","primary vertex y [cm]","N evts",200,-0.21,0.19);
  TH1D *hPVz = MakeHisto("hPVz","primary vertex z","primary vertex z [cm]","N evts",200,-20.05,19.95);
  TH1D *hPVrho = MakeHisto("hPVrho","primary vertex rho","primary vertex #rho [cm]","N evts",2000,0,1);
  TH1D *hPVchi2 = MakeHisto("hPVchi2","primary vertex chi2","primary vertex #chi^{2}","N evts",500,0,100);
  TH1D *hPVndof = MakeHisto("hPVndof","primary vertex ndof","primary vertex ndof","N evts",100,0,100);
  TH1D *hPVchi2ndof = MakeHisto("hPVchi2ndof","primary vertex chi2/ndof","primary vertex #chi^{2}/ndof","N evts",200,0,10); 
  TH1D *hPVnTrack = MakeHisto("hPVnTrack","n track primary vertex","n track primary vertex","N evts",100,0,100);
  TH1D *hPVfTrack = MakeHisto("hPVfTrack","fraction of tracks primary vertex","fraction of tracks primary vertex","N evts",100,0,100);
  TH1D *hPVTrackPt = MakeHisto("hPVTrackPt","primary vertex track pT","track pT","N tracks",150,0,15);
  TH1D *hPVTrackEta = MakeHisto("hPVTrackEta","primary vertex track eta","track #eta","N tracks",120,-3,3);
  
  cout<<endl;

  //-- second vertices
  TH1D *hSVx = MakeHisto("hSVx","second vertex x","second vertex x [cm]","N evts",200,-0.21,0.19);
  TH1D *hSVy = MakeHisto("hSVy","second vertex y","second vertex y [cm]","N evts",200,-0.21,0.19);
  TH1D *hSVz = MakeHisto("hSVz","second vertex z","second vertex z [cm]","N evts",200,-20.05,19.95);
  TH1D *hSVrho = MakeHisto("hSVrho","second vertex rho","second vertex #rho [cm]","N evts",2000,0,1);
  TH1D *hSVchi2 = MakeHisto("hSVchi2","second vertex chi2","second vertex #chi^{2}","N evts",500,0,100);
  TH1D *hSVndof = MakeHisto("hSVndof","second vertex ndof","second vertex ndof","N evts",100,0,100);
  TH1D *hSVchi2ndof = MakeHisto("hSVchi2ndof","second vertex chi2/ndof","second vertex #chi^{2}/ndof","N evts",200,0,10); 
  TH1D *hSVnTrack = MakeHisto("hSVnTrack","n track second vertex","n track second vertex","N evts",100,0,100);
  TH1D *hSVfTrack = MakeHisto("hSVfTrack","fraction of tracks second vertex","fraction of tracks second vertex","N evts",100,0,100);
  cout<<endl;

  //-- vertices correlation
  TH1D *hDeltaPVzSVz = MakeHisto("hDeltaPVzSVz","PV z - SV z","delta z between primary and second vertex","N evts",1200,-30,30);
  TH2D *hPVzSVz = MakeHisto("hPVzSVz","SV z versus PV z","primary vertex z [cm]","second vertex z [cm]",200,-20.05,19.95,200,-20.05,19.95);
  TH2D *hPVnTrackSVnTrack = MakeHisto("hPVnTrackSVnTrack","SV n track versus PV n track","n track primary vertex","n track second vertex",
				      100,0,100,100,0,100);
  TH2D *hPVfTrackSVfTrack = MakeHisto("hPVfTrackSVfTrack","SV f track versus PV f track","fraction of tracks primary vertex",
				      "fraction of tracks second vertex",100,0,100,100,0,100);

  TH2D *hSVnTrackDeltaPVzSVz = MakeHisto("hSVnTrackDeltaPVzSVz","SV n track versus delta z between primary and second vertex",
					 "delta z between primary and second vertex","SV n track",300,-30,30,120,0,120);
  cout<<endl;

  //-- pileup information
  TH1D *hnbxint = MakeHisto("hnbxint","in time bunch crossing multiplicity","in time bunch crossing multiplicity","number of events",11,-0.5,10.5);
  TH1D *hnbxearly = MakeHisto("hnbxearly","early bunch crossing multiplicity","early bunch crossing multiplicity","number of events",11,-0.5,10.5);
  TH1D *hnbxlate = MakeHisto("hnbxlate","late bunch crossing multiplicity","late bunch crossing multiplicity","number of events",11,-0.5,10.5);
  TH1D *hnbxtot = MakeHisto("hnbxtot","total bunch crossing multiplicity","total bunch crossing multiplicity","number of events",11,-0.5,10.5);  

  TH1D *hPUint = MakeHisto("hPUint","1 + in time PU multiplicity","1 + in time PU multiplicity","number of events",11,-0.5,10.5);
  TH1D *hPUearly = MakeHisto("hPUearly","1 + early PU multiplicity","1 + early PU multiplicity","number of events",11,-0.5,10.5);
  TH1D *hPUlate = MakeHisto("hPUlate","1 + late PU multiplicity","1 + late PU multiplicity","number of events",11,-0.5,10.5);
  TH1D *hPUtot = MakeHisto("hPUtot","1 + total PU multiplicity","1 + total PU multiplicity","number of events",11,-0.5,10.5);

  TH1D* hNumInt = MakeHisto("hNumInt","interactions multiplicity","interactions multiplicity","number of events",100,0,10); 

  TH2D* hPUintVNum = MakeHisto("hPUintVNum","vertex multiplicity versus 1 + in time PU multiplicity","1 + in time PU multiplicity","vertex multiplicity",
			       11,-0.5,10.5,11,-0.5,10.5);
  TH2D* hPUearlyVNum = MakeHisto("hPUearlyVNum","vertex multiplicity versus 1 + early PU multiplicity","1 + early PU multiplicity","vertex multiplicity",
				 11,-0.5,10.5,11,-0.5,10.5);
  TH2D* hPUlateVNum = MakeHisto("hPUlateVNum","vertex multiplicity versus 1 + late PU multiplicity","1 + late PU multiplicity","vertex multiplicity",
				11,-0.5,10.5,11,-0.5,10.5);
  TH2D* hPUtotVNum = MakeHisto("hPUtotVNum","vertex multiplicity versus 1 + total PU multiplicity","1 + total PU multiplicity","vertex multiplicity",
			       11,-0.5,10.5,11,-0.5,10.5);      
  TH2D* hPUintVNum0 = MakeHisto("hPUintVNum0","vertex multiplicity versus 1 + in time PU multiplicity","1 + in time PU multiplicity","vertex multiplicity",
				11,-0.5,10.5,11,-0.5,10.5);

  //-- weight  	    
  double weight = 1;

  //-- vertex weight data
  double weight_vertex = 1;

  //-- PU weight data
  double weight_pu = 1;

  //-- PU weight monte carlo
  double weight_pu_array[10];
  for(int i = 0; i < 10; i++) weight_pu_array[i] = 1;
  int Niteration_pu = 4;
  if(!isData && pu_reweight) GetPUWeight(AllMC,weight_pu_array,Niteration_pu);

  //-- trigger 
  std::vector<std::string> HLTrequested;

  for(int i = 0; i < 8; ++i) {
    string temp = "HLT_ZeroBias_part" + std::to_string(i) +"_v1";
    HLTrequested.push_back(temp);
  }

  std::vector<std::string> L1TTrequested;
  L1TTrequested.push_back("L1Tech_BPTX_plus_AND_minus.v0");

  //-- Determine total number of files and events 
  TIter temp_next(filelist); 
  TObjString* temp_itfile = 0;
  
  int Nfile_all = 0;
  int Nevt_all_files = 0;
 
  int file_nb_max = -1;  
  int file_nb_read = 0;

  while((temp_itfile = (TObjString*)temp_next()) && (file_nb_read < file_nb_max || file_nb_max == -1)) {
    
    file_nb_read++;
    
    TFile* file = TFile::Open(inputdir+temp_itfile->GetString(),"READ");
   
    if (!file) {
      cout<<"Error in TrackAnalyzer: could not open file "<<temp_itfile->GetString()<<endl;
      continue;
    }    
    
    TTree *tree = new TTree("TrackAnalysis","");
    file->GetObject("analysis/TrackAnalysis",tree);
   
    Nfile_all++;
    Nevt_all_files+=tree->GetEntriesFast();
    delete tree;
    file->Close();
    delete file;
  }
  
  cout<<endl<<"you are going to loop on "<<Nfile_all<<" files containing "<<Nevt_all_files<<" events"<<endl<<endl;
  cout<<"press enter to continue"<<endl;
  getchar();

  //-- Loop over the different files
  int file_nb = 0;	     

  TIter next(filelist);
  TObjString* itfile = 0;
  TString filename = "";

  while((itfile = (TObjString*)next()) && (file_nb < file_nb_max || file_nb_max == -1)) {
    
    file_nb++;
    
    filename.Clear();
    filename = itfile->GetString();
    
    cout<<endl<<"open file "<<file_nb<<" with name: "<<itfile->GetString()<<endl;
    
    TFile* file = TFile::Open(inputdir+itfile->GetString(),"READ");
    
    if (!file) { 
      cout<<"Error in TrackAnalyzer: could not open file "<<itfile->GetString()<<endl; 
      continue; 
    } 
    
    //-- Get tree from file
    TTree *tree = new TTree("TrackAnalysis","");
    file->GetObject("analysis/TrackAnalysis",tree);
    
    //-- Get Branches
    GetBranches(tree,isData);
    
    //-- Set Branch Addresses
    SetBranchAddresses(tree,isData);
    
    int Nevents = tree->GetEntriesFast();
    cout <<"number of events in file "<<file_nb<<" = "<<Nevents<<endl<<endl;
    
    //-- Loop over the events 
    for (int ievt=0;ievt<Nevents;ievt++) {

      //-- Get Entries
      GetTreeEntries(ievt,isData);

      HFplusTowerLeadingEnergy*=HFscaling;
      HFminusTowerLeadingEnergy*=HFscaling; 

      //-- rho with respect to the beam spot
      double rhoBS[nVertex];
      double VxBS[nVertex];
      double VyBS[nVertex];

      for (int ivert = 0; ivert < nVertex; ++ivert) {	
	double deltax = Vx->at(ivert) - BSx;
	double deltay = Vy->at(ivert) - BSy;
	VxBS[ivert] = deltax;
	VyBS[ivert] = deltay;
	rhoBS[ivert] = TMath::Sqrt(deltax*deltax+deltay*deltay);
      }

      if(nVertex > 1 && BScheck) {
	for (int ivert = 0; ivert < nVertex; ++ivert) {	
	  cout<<"vertex "<<ivert+1<<endl<<endl;
	  cout<<"Vx = "<<Vx->at(ivert)<<endl;
	  cout<<"Vy = "<<Vy->at(ivert)<<endl<<endl;
	  cout<<"BSx = "<<BSx<<endl;
	  cout<<"BSy = "<<BSy<<endl<<endl;
	  cout<<"rho = "<<Vrho->at(ivert)<<endl;
	  cout<<"rho BS = "<<rhoBS[ivert]<<endl<<endl;
	  getchar();
	}
      }
      
      //-- vertex information - at least one good vertex
      bool at_least_one_gv = false;
      bool at_least_one_gv_20 = false;

      bool good_vertex[int(nVertex)];
      for (int ivertex = 0; ivertex < nVertex; ++ivertex) good_vertex[ivertex] = false;

      int ngood_vertex = 0;

      for (int ivertex = 0; ivertex < nVertex; ++ivertex) {	
	
	if(VisFake->at(ivertex)) continue;
	if(!VisValid->at(ivertex)) continue;
	
	if(rhoBS[ivertex] > 0.2) continue;  //-- rhoBS <= 0.2

	int ngood_track = 0;
	
	for (int itrack = 0; itrack < VnTrack->at(ivertex); itrack++) {  
	  
	  if(VTrackQuality->at(ivertex).at(itrack) != 2) continue; //-- high quality
		  
	  double pTerror = 1000;
	  if(VTrackPt->at(ivertex).at(itrack) != 0) 
	    pTerror = 100*VTrackPtError->at(ivertex).at(itrack)/VTrackPt->at(ivertex).at(itrack);
	  if(pTerror > 10) continue; //-- pTerror <= 10
	
	  double dxysigmaxy = 1000;
	  if(VTrackdxyError->at(ivertex).at(itrack) != 0)
	    dxysigmaxy = VTrackdxy->at(ivertex).at(itrack)/VTrackdxyError->at(ivertex).at(itrack);
	  if(dxysigmaxy > 3) continue; //-- dxysigmaxy <= 3
	
	  double dzsigmaz = 1000;	  
	  if(VTrackdzError->at(ivertex).at(itrack) != 0) 
	    dzsigmaz = VTrackdz->at(ivertex).at(itrack)/VTrackdzError->at(ivertex).at(itrack);
	  if(dzsigmaz > 3) continue; //-- dzsigmaz <= 3

	  int npixelhit = VTrackValidPixelHit->at(ivertex).at(itrack);
	  double abseta = TMath::Abs(VTrackEta->at(ivertex).at(itrack));

	  if(npixelhit < 2) continue;  //-- n pixel hits >= 2
	  if(npixelhit < 3 && pixel_cut) continue;  //-- n pixel hits >= 3

	  if(npixelhit < 3 && abseta < 1 && combined_cut) continue;
	  if(npixelhit < 2 && abseta > 1 && combined_cut) continue;
	  ngood_track++;
	} //-- end loop over associated tracks
	
	if(ngood_track < 2) continue; //-- nTrack >= 2

	if(TMath::Abs(Vz->at(ivertex)) > 20) continue;   //-- |Vz| <= 20                                                       
	at_least_one_gv_20 = true;

	if(TMath::Abs(Vz->at(ivertex)) > 15) continue;   //-- |Vz| <= 15
	at_least_one_gv = true;
	good_vertex[ivertex] = true;
	ngood_vertex++;
      } //-- end loop over vertex
      
      //-- vertex information - exactly one good vertex
      bool exactly_one_gv = false;					      
      if(at_least_one_gv == true && nVertex == 1) exactly_one_gv = true;

      bool exactly_one_gv_20 = false;
      if(at_least_one_gv_20 == true && nVertex == 1) exactly_one_gv_20 = true;

      //-- vertex information - good vertex collection
      double Vx_gv[int(ngood_vertex)];
      double Vy_gv[int(ngood_vertex)];
      double Vz_gv[int(ngood_vertex)];
      double Vrho_gv[int(ngood_vertex)];     

      double Vchi2_gv[int(ngood_vertex)];	  
      double Vndof_gv[int(ngood_vertex)];	   
      double Vchi2ndof_gv[int(ngood_vertex)];	  
      
      double Vntrack_gv[int(ngood_vertex)]; 
      for (int igood = 0; igood < ngood_vertex; ++igood) Vntrack_gv[igood] = 0;
      double ntrack_all_gv = 0; 

      int igood = 0;
      for (int ivertex = 0; ivertex < nVertex; ++ivertex) {	

	if(good_vertex[ivertex] == false) continue; //-- only consider good vertices
	
	Vx_gv[igood] = Vx->at(ivertex);
	Vy_gv[igood] = Vy->at(ivertex);
	Vz_gv[igood] = Vz->at(ivertex);
	Vrho_gv[igood] = Vrho->at(ivertex);

	Vchi2_gv[igood] = Vchi2->at(ivertex);	       
	Vndof_gv[igood] = Vndof->at(ivertex);	        
	Vchi2ndof_gv[igood] = Vchi2ndof->at(ivertex);
	
	//-- determine number of high quality tracks with pT > 0.5 GeV attached to each vertex 
	for (int itrack = 0; itrack < VnTrack->at(ivertex); itrack++) {  	    
	  if(VTrackPt->at(ivertex).at(itrack) < 0.5) continue;       //-- cut on track pT
	  if(VTrackQuality->at(ivertex).at(itrack) != 2) continue;   //-- cut on track quality

	  int npixelhit = VTrackValidPixelHit->at(ivertex).at(itrack);
	  double abseta = TMath::Abs(VTrackEta->at(ivertex).at(itrack));

	  if(npixelhit < 2) continue;  //-- cut on pixel hit
	  if(npixelhit < 3 && pixel_cut) continue;  //-- cut on pixel hit

	  if(npixelhit < 3 && abseta < 1 && combined_cut) continue;
	  if(npixelhit < 2 && abseta > 1 && combined_cut) continue;
	  Vntrack_gv[igood]++;
	}
         
	ntrack_all_gv+=Vntrack_gv[igood];
	igood++;
      }
	 
      if(GoodVertex) {
	cout<<"number of vertices = "<<nVertex<<endl;
	cout<<"number of good vertices = "<<ngood_vertex<<endl<<endl;
	
	int igv = 0;

	for (int ivertex = 0; ivertex < nVertex; ++ivertex) {
	  if(good_vertex[ivertex] == true) cout<<"vertex "<<ivertex+1<<": good"<<endl;
	  if(good_vertex[ivertex] == false) cout<<"vertex "<<ivertex+1<<": NOT good"<<endl;

	  if(good_vertex[ivertex] == false){
	    cout<<"Vx = "<<Vx->at(ivertex)<<endl;
	    cout<<"Vy = "<<Vy->at(ivertex)<<endl;
	    cout<<"Vz = "<<Vz->at(ivertex)<<endl;
	    cout<<"Vrho = "<<Vrho->at(ivertex)<<endl;
	    cout<<"Vchi2 = "<<Vchi2->at(ivertex)<<endl;	       
	    cout<<"Vndof = "<<Vndof->at(ivertex)<<endl;	        
	    cout<<"Vchi2ndof = "<<Vchi2ndof->at(ivertex)<<endl<<endl;
	  }

	  if(good_vertex[ivertex] == true){
	    cout<<"Vx = "<<Vx->at(ivertex)<<" - Vx good = "<<Vx_gv[igv]<<endl;		    
	    cout<<"Vy = "<<Vy->at(ivertex)<<" - Vy good = "<<Vy_gv[igv]<<endl;		    
	    cout<<"Vz = "<<Vz->at(ivertex)<<" - Vz good = "<<Vz_gv[igv]<<endl;		    
	    cout<<"Vrho = "<<Vrho->at(ivertex)<<" - Vrho good = "<<Vrho_gv[igv]<<endl;	    
	    cout<<"Vchi2 = "<<Vchi2->at(ivertex)<<" - Vchi2 good = "<<Vchi2_gv[igv]<<endl;	    
	    cout<<"Vndof = "<<Vndof->at(ivertex)<<" - Vndof good = "<<Vndof_gv[igv]<<endl;	    
	    cout<<"Vchi2ndof = "<<Vchi2ndof->at(ivertex)<<" - Vchi2ndof good = "<<Vchi2ndof_gv[igv]<<endl;
	    cout<<"Vntrack good = "<<Vntrack_gv[igv]<<endl<<endl;
	    igv++;
	  }
	}

	cout<<endl;
	getchar();

      } //-- GoodVertex check

      //-- vertex weight montecarlo 
      weight_vertex = 1; 
      int Niteration_vertex = 2;
      if(!isData && exactly_one_gv_20 && vz_reweight) weight_vertex = GetVertexWeight(Vz->at(0),AllMC,Niteration_vertex); //-- otherwise set to 1

      if(VertexDebug) {                       
	cout<<"n vertex = "<<nVertex<<" - vertex weight = "<<weight_vertex<<endl;
	getchar();
      }

      //-- pu weight monte carlo
      weight_pu = 1;
      if(!isData && pu_reweight) weight_pu = weight_pu_array[int(PUint)];

      //-- global weight
      weight = weight_pu*weight_vertex;

      if(PUdebug) {
	cout<<"in time PU multiplicity = "<<PUint<<" - weight PU = "<<weight_pu<<endl;
	cout<<"weight = "<<weight<<endl;
	getchar();
      }

      //-- vertex reconstruction efficiency
      hPUintVNum0->Fill(1+PUint,ngood_vertex,weight);

      //-- number of events 
      Nevt_tot+=weight;
      Nevt_tot_no_weight+=1;
      Nevt_tot_no_pu_weight+=weight_vertex;
      Nevt_tot_no_vz_weight+=weight_pu; 

      if(isData) hselection_data->Fill(1,weight);
      if(!isData) hselection_moca_reco->Fill(1,weight);

      if((int (Nevt_tot_no_weight))%20000 == 0) {
	cout<<"total number of events done - no weight = "<<Nevt_tot_no_weight<<" ("<<100*Nevt_tot_no_weight/Nevt_all_files<<"%)"<<endl;
	cout<<"total number of events done - global weight = "<<Nevt_tot<<endl<<endl;
	if(isData) cout<<"number of selected events in data - global weight = "<<Ndata_sel<<endl<<endl;
	if(!isData) cout<<"number of selected events in mc at reco level - global weight = "<<Nmoca_reco_sel<<endl<<endl;
      }

      //-- generated level information (no pT cut, no reweighting)								  

      if(!isData) {

	double ptgen_max = 0;
	double etagen_max = 0;
	double phigen_max = 0;

	int numch = 0;          
	int numch500 = 0;
        int numhad = 0;
  
	int ngenparticle = ParticleEnergy->size();

	//-- loop over generated particles
	for(int igen = 0; igen < ngenparticle; ++igen) {
	  
	  if(ParticleStatus->at(igen) != 1) continue;
	  if(ParticleCharge->at(igen) == 0) continue;
	  if(TMath::Abs(ParticleEta->at(igen)) > 2.4) continue;
	  numch++;

	  hChPt->Fill(ParticlePt->at(igen));
	  hChEta->Fill(ParticleEta->at(igen));
	  hChPhi->Fill(ParticlePhi->at(igen));

	  int id = ParticleId->at(igen);
	  if(isProton(AllMC,id) || isPion(AllMC,id) || isKaon(AllMC,id)) {
	    numhad++;
	    hHadPt->Fill(ParticlePt->at(igen));  
	    hHadEta->Fill(ParticleEta->at(igen));
	    hHadPhi->Fill(ParticlePhi->at(igen));	
	  }  
	  
	  if(ParticlePt->at(igen) > ptgen_max) {
	    ptgen_max = ParticlePt->at(igen);
	    etagen_max = ParticleEta->at(igen);
	    phigen_max = ParticlePhi->at(igen);
	  }

	  if(ParticlePt->at(igen) < 0.5) continue;
	  numch500++;

	  h500ChPt->Fill(ParticlePt->at(igen));  
	  h500ChEta->Fill(ParticleEta->at(igen));
	  h500ChPhi->Fill(ParticlePhi->at(igen));	  


	} //-- end loop generated particles

	hChNum->Fill(numch);
	h500ChNum->Fill(numch500);
	hHadNum->Fill(numhad);

	if(ptgen_max > 0) {
	  hChLeadingPt->Fill(ptgen_max); 
	  hChLeadingEta->Fill(etagen_max);
	  hChLeadingPhi->Fill(phigen_max);
	}                                   

	if(ptgen_max > 0.5) {		    
	  h500ChLeadingPt->Fill(ptgen_max);    
	  h500ChLeadingEta->Fill(etagen_max);  
	  h500ChLeadingPhi->Fill(phigen_max);  
	}                                   
      
      } //-- end generated level information (no pT cut, no reweighting)									       
       
			 
      //-- Filter on data       
      bool run_selection = false;
      bool L1TT_selection = false;
      bool HLT_selection = false;
      bool data_sel_at_least_one_gv = false;
      bool data_sel = false;
      
      if(isData) {
	
	//-- filter on run and lumi section 	
	if(Run == 251721 && LumiSection >= 90) run_selection = true;
	
	//-- filter on L1TT trigger
	for(unsigned int il1tt = 0; il1tt < L1TTname->size(); ++il1tt) {
	  for(unsigned int irequest = 0; irequest < L1TTrequested.size(); ++irequest) {
	    if(L1TTname->at(il1tt).compare(L1TTrequested.at(irequest)) == 0) {
	      if(L1TTdecision->at(il1tt)) {
		//-- cout<<"requested L1TT path "<<L1TTname->at(il1tt)<<" is present with decision "<<L1TTdecision->at(il1tt)<<endl;
		L1TT_selection = true;
	      }
	    }
	  }
	}	      
	
	//-- filter on HLT trigger
	for(unsigned int ihlt = 0; ihlt < HLTname->size(); ++ihlt) {
	  for(unsigned int irequest = 0; irequest < HLTrequested.size(); ++irequest) {
	    if(HLTname->at(ihlt).compare(HLTrequested.at(irequest)) == 0) {
	      if(HLTdecision->at(ihlt)) {
		//-- cout<<"requested HLT path "<<HLTname->at(ihlt)<<" is present with decision "<<HLTdecision->at(ihlt)<<endl;
		HLT_selection = true;
	      }
	    }
	  }
	}
	
	if(run_selection) hselection_data->Fill(2,weight);
	if(run_selection && L1TT_selection) hselection_data->Fill(3,weight); 
	if(run_selection && L1TT_selection && HLT_selection) hselection_data->Fill(4,weight);

	if(run_selection && L1TT_selection && HLT_selection && at_least_one_gv) {
	  hselection_data->Fill(5,weight);
	  data_sel_at_least_one_gv = true;
	}

	if(run_selection && L1TT_selection && HLT_selection && at_least_one_gv && exactly_one_gv) {
	  hselection_data->Fill(6,weight);
	  data_sel = true;

	  Ndata_sel+=weight;									   
	  Ndata_sel_no_weight+=1;
	  Ndata_sel_no_pu_weight+=weight_vertex;
	  Ndata_sel_no_vz_weight+=weight_pu;
	}
	
      } //-- end loop on data
      
      //-- Filter on moca reco  
      bool moca_reco_sel = false;
      bool moca_reco_sel_at_least_one_gv = false;
      
      if(!isData) {
	if(at_least_one_gv) {
	  hselection_moca_reco->Fill(2,weight);
	  moca_reco_sel_at_least_one_gv = true;
	}
	
	if(at_least_one_gv && exactly_one_gv) {
	  hselection_moca_reco->Fill(3,weight);
	  moca_reco_sel = true;

	  Nmoca_reco_sel+=weight;									  
	  Nmoca_reco_sel_no_weight+=1;
	  Nmoca_reco_sel_no_pu_weight+=weight_vertex;
	  Nmoca_reco_sel_no_vz_weight+=weight_pu;
	}
      } //-- end loop on moca reco

      //-- vertices and tracks - at least one good vertex
                                                                                                                            
      double fTrackPV = 0;
      double fTrackSV = 0;      

      vector <double> ChPt;			          
      vector <double> ChEta;			    
      vector <double> ChPhi;			    
      
      vector <double> MChPt;			    
      vector <double> MChEta;			    
      vector <double> MChPhi;			    
      
      vector <double> MTPt; 			    
      vector <double> MTEta;			    
      vector <double> MTPhi;       		    
      
      int nch_matching = 0;				    
      int ntrack_matching = 0;                            

      //-- Apply filter data - moca reco - at least one good vertex
      if((isData && data_sel_at_least_one_gv) || (!isData && moca_reco_sel_at_least_one_gv)) {
	
	if(isData && LumiSection < 90) cout<<"warning: "<<"LumiSection = "<<LumiSection<<" < 90"<<endl;
	
	//-- pileup information												     
	if(!isData) {
	  hnbxint->Fill(nbxint,weight); 
	  hnbxearly->Fill(nbxearly,weight); 
	  hnbxlate->Fill(nbxlate,weight); 
	  hnbxtot->Fill(nbxtot,weight); 

	  hPUint->Fill(1+PUint,weight); 
	  hPUearly->Fill(1+PUearly,weight); 
	  hPUlate->Fill(1+PUlate,weight); 
	  hPUtot->Fill(1+PUtot,weight); 

	  hNumInt->Fill(NumInteraction,weight);

	  hPUintVNum->Fill(1+PUint,ngood_vertex,weight);
	  hPUearlyVNum->Fill(1+PUearly,ngood_vertex,weight);
	  hPUlateVNum->Fill(1+PUlate,ngood_vertex,weight);
	  hPUtotVNum->Fill(1+PUtot,ngood_vertex,weight);     
	}                                           

	hVNum->Fill(nVertex,weight);
	hVNumGood->Fill(ngood_vertex,weight);
	hVNum2D->Fill(nVertex,ngood_vertex,weight);

	if(ngood_vertex >= 2) {

      	  //-- primary vertices
	  hPVx->Fill(Vx_gv[0],weight);
	  hPVy->Fill(Vy_gv[0],weight);
	  hPVz->Fill(Vz_gv[0],weight);
	  hPVrho->Fill(Vrho_gv[0],weight);
	  hPVchi2->Fill(Vchi2_gv[0],weight);
	  hPVndof->Fill(Vndof_gv[0],weight);
	  hPVchi2ndof->Fill(Vchi2ndof_gv[0],weight);
	  hPVnTrack->Fill(Vntrack_gv[0],weight);   

	  for (int itrack = 0; itrack < VnTrack->at(0); itrack++) {	      

	    if(VTrackPt->at(0).at(itrack) < 0.5) continue; //-- cut on track pT
	    if(VTrackQuality->at(0).at(itrack) != 2) continue; //-- cut on track quality

	    int npixelhit = VTrackValidPixelHit->at(0).at(itrack);
	    double abseta = TMath::Abs(VTrackEta->at(0).at(itrack));

	    if(npixelhit < 2) continue;  //-- cut on pixel hit
	    if(npixelhit < 3 && pixel_cut) continue;  //-- cut on pixel hit
	    
	    if(npixelhit < 3 && abseta < 1 && combined_cut) continue;
	    if(npixelhit < 2 && abseta > 1 && combined_cut) continue;

	    hPVTrackPt->Fill(VTrackPt->at(0).at(itrack),weight);
	    hPVTrackEta->Fill(VTrackEta->at(0).at(itrack),weight);
	  }
	  
	  //-- second vertices
	  hSVx->Fill(Vx_gv[1],weight);
	  hSVy->Fill(Vy_gv[1],weight);
	  hSVz->Fill(Vz_gv[1],weight);
	  hSVrho->Fill(Vrho_gv[1],weight);
	  hSVchi2->Fill(Vchi2_gv[1],weight);
	  hSVndof->Fill(Vndof_gv[1],weight);
	  hSVchi2ndof->Fill(Vchi2ndof_gv[1],weight);
	  hSVnTrack->Fill(Vntrack_gv[1],weight);
	    
	  hPVzSVz->Fill(Vz_gv[0],Vz_gv[1],weight);
	  hDeltaPVzSVz->Fill(Vz_gv[0]-Vz_gv[1],weight);
	  hSVnTrackDeltaPVzSVz->Fill(Vz_gv[0]-Vz_gv[1],Vntrack_gv[1],weight); 
	  hPVnTrackSVnTrack->Fill(Vntrack_gv[0],Vntrack_gv[1],weight);
 
	  //-- fraction of tracks
	  if(ntrack_all_gv != 0) {
	    fTrackPV = 100*Vntrack_gv[0]/ntrack_all_gv;
	    hPVfTrack->Fill(fTrackPV,weight);

	    fTrackSV = 100*Vntrack_gv[1]/ntrack_all_gv;
	    hSVfTrack->Fill(fTrackSV,weight);

	    hPVfTrackSVfTrack->Fill(fTrackPV,fTrackSV,weight);
	  } //-- end fraction of tracks
	  
	} //-- end ngood_vertex >= 2
      } //-- end loop on data - moca reco at least one good vertex

      //-- Apply filter data - moca reco - exactly one good vertex
      if((isData && data_sel) || (!isData && moca_reco_sel)) {
	
	if(isData && LumiSection < 90) cout<<"warning: "<<"LumiSection = "<<LumiSection<<" < 90"<<endl;

	//-- HF tower 
	int nHFTower_HFcut[nHFcut];
	for(int ihf = 0; ihf < nHFcut; ++ihf) nHFTower_HFcut[ihf] = 0;
	
	double HFTowerEtot_HFcut[nHFcut];
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFTowerEtot_HFcut[ihf] = 0;
	
	double HFTowerLeadingEnergy_HFcut[nHFcut];
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFTowerLeadingEnergy_HFcut[ihf] = 0;

	//-- HF plus tower 								  
	int nHFplusTower_HFcut[nHFcut];						  
	for(int ihf = 0; ihf < nHFcut; ++ihf) nHFplusTower_HFcut[ihf] = 0;		  
	
	double HFplusTowerEtot_HFcut[nHFcut];					  
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFplusTowerEtot_HFcut[ihf] = 0;	  
	
	double HFplusTowerLeadingEnergy_HFcut[nHFcut];				  
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFplusTowerLeadingEnergy_HFcut[ihf] = 0;	

	//-- HF minus tower 								
	int nHFminusTower_HFcut[nHFcut];						  	
	for(int ihf = 0; ihf < nHFcut; ++ihf) nHFminusTower_HFcut[ihf] = 0;		
	
	double HFminusTowerEtot_HFcut[nHFcut];					  	
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFminusTowerEtot_HFcut[ihf] = 0;	  	
	
	double HFminusTowerLeadingEnergy_HFcut[nHFcut];				  	
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFminusTowerLeadingEnergy_HFcut[ihf] = 0;	

	//-- no cut
	int nHFplusTowerNoCut = 0;						  
	int nHFminusTowerNoCut = 0;						  	

	for(int itower = 0; itower < nHFTower; ++itower) {

	  if(HFTowerEta->at(itower) > 0) { 		
	    hHFplusTowerEnergyNoCut->Fill(HFscaling*HFTowerEnergy->at(itower),weight);
	    nHFplusTowerNoCut++;						  
	  }

	  if(HFTowerEta->at(itower) < 0) {
	    hHFminusTowerEnergyNoCut->Fill(HFscaling*HFTowerEnergy->at(itower),weight);
	    nHFminusTowerNoCut++;						  	
	  }

	}

	//cout<<"nHFplusTowerNoCut = "<<nHFplusTowerNoCut<<endl;
	//cout<<"nHFminusTowerNoCut = "<<nHFminusTowerNoCut<<endl;

	hnHFplusTowerNoCut->Fill(nHFplusTowerNoCut,weight);										   
	hnHFminusTowerNoCut->Fill(nHFminusTowerNoCut,weight);										   

	//-- different cuts
	for(int ihf = 0; ihf < nHFcut; ++ihf) {
	  
	  for(int itower = 0; itower < nHFTower; ++itower) {

	    if(HFscaling*HFTowerEnergy->at(itower) > EHFcut[ihf]) {
	      
	      //-- HF plus and HF minus together
	      hHFTowerEnergy[ihf]->Fill(HFscaling*HFTowerEnergy->at(itower),weight);
	      hHFTowerEta[ihf]->Fill(HFTowerEta->at(itower),weight);
	      hHFTowerPhi[ihf]->Fill(HFTowerPhi->at(itower),weight);
	      hHFTowerEtaPhi[ihf]->Fill(HFTowerEta->at(itower),HFTowerPhi->at(itower),weight);

	      nHFTower_HFcut[ihf]++;
	      HFTowerEtot_HFcut[ihf]+=HFscaling*HFTowerEnergy->at(itower);
	      
	      if(HFscaling*HFTowerEnergy->at(itower) >  HFTowerLeadingEnergy_HFcut[ihf]) 
		HFTowerLeadingEnergy_HFcut[ihf] = HFscaling*HFTowerEnergy->at(itower);

	      //-- HF plus alone
	      if(HFTowerEta->at(itower) > 0) {		
		hHFplusTowerEnergy[ihf]->Fill(HFscaling*HFTowerEnergy->at(itower),weight);
		hHFplusTowerEta[ihf]->Fill(HFTowerEta->at(itower),weight);
		hHFplusTowerPhi[ihf]->Fill(HFTowerPhi->at(itower),weight);
		hHFplusTowerEtaPhi[ihf]->Fill(HFTowerEta->at(itower),HFTowerPhi->at(itower),weight);

		nHFplusTower_HFcut[ihf]++;
		HFplusTowerEtot_HFcut[ihf]+=HFscaling*HFTowerEnergy->at(itower);
	      
		if(HFscaling*HFTowerEnergy->at(itower) >  HFplusTowerLeadingEnergy_HFcut[ihf]) 
		  HFplusTowerLeadingEnergy_HFcut[ihf] = HFscaling*HFTowerEnergy->at(itower);     
	      }

	      //-- HF minus alone
	      if(HFTowerEta->at(itower) < 0) {		
		hHFminusTowerEnergy[ihf]->Fill(HFscaling*HFTowerEnergy->at(itower),weight);
		hHFminusTowerEta[ihf]->Fill(HFTowerEta->at(itower),weight);
		hHFminusTowerPhi[ihf]->Fill(HFTowerPhi->at(itower),weight);
		hHFminusTowerEtaPhi[ihf]->Fill(HFTowerEta->at(itower),HFTowerPhi->at(itower),weight);

		nHFminusTower_HFcut[ihf]++;
		HFminusTowerEtot_HFcut[ihf]+=HFscaling*HFTowerEnergy->at(itower);
		
		if(HFscaling*HFTowerEnergy->at(itower) >  HFminusTowerLeadingEnergy_HFcut[ihf]) 
		  HFminusTowerLeadingEnergy_HFcut[ihf] = HFscaling*HFTowerEnergy->at(itower);
	      }
 
	    } //-- end EHFcut
	  } //-- end loop over HFtower

	  //-- HF plus and HF minus together
	  hnHFTower[ihf]->Fill(nHFTower_HFcut[ihf],weight);
	  hHFTowerEtot[ihf]->Fill(HFTowerEtot_HFcut[ihf],weight);
	  if(HFTowerLeadingEnergy_HFcut[ihf] > 0) hHFTowerLeadingEnergy[ihf]->Fill(HFTowerLeadingEnergy_HFcut[ihf],weight);

	  //-- HF plus alone
	  hnHFplusTower[ihf]->Fill(nHFplusTower_HFcut[ihf],weight);			   
	  hHFplusTowerEtot[ihf]->Fill(HFplusTowerEtot_HFcut[ihf],weight);		   
	  if(HFplusTowerLeadingEnergy_HFcut[ihf] > 0) hHFplusTowerLeadingEnergy[ihf]->Fill(HFplusTowerLeadingEnergy_HFcut[ihf],weight);

	  //-- HF minus alone								   
	  hnHFminusTower[ihf]->Fill(nHFminusTower_HFcut[ihf],weight);			   
	  hHFminusTowerEtot[ihf]->Fill(HFminusTowerEtot_HFcut[ihf],weight);		   
	  if(HFminusTowerLeadingEnergy_HFcut[ihf] > 0) hHFminusTowerLeadingEnergy[ihf]->Fill(HFminusTowerLeadingEnergy_HFcut[ihf],weight);
	
	} //-- end loop over nHFcut

	//-- HF plus rechit 								  
	int nHFplusRechit_HFcut[nHFcut];						  
	for(int ihf = 0; ihf < nHFcut; ++ihf) nHFplusRechit_HFcut[ihf] = 0;		  
	
	double HFplusRechitEtot_HFcut[nHFcut];					  
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFplusRechitEtot_HFcut[ihf] = 0;	  
	
	double HFplusRechitLeadingEnergy_HFcut[nHFcut];				  
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFplusRechitLeadingEnergy_HFcut[ihf] = 0;	

	//-- HF minus rechit 								  
	int nHFminusRechit_HFcut[nHFcut];						  
	for(int ihf = 0; ihf < nHFcut; ++ihf) nHFminusRechit_HFcut[ihf] = 0;		  
	
	double HFminusRechitEtot_HFcut[nHFcut];					  
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFminusRechitEtot_HFcut[ihf] = 0;	  
	
	double HFminusRechitLeadingEnergy_HFcut[nHFcut];				  
	for(int ihf = 0; ihf < nHFcut; ++ihf) HFminusRechitLeadingEnergy_HFcut[ihf] = 0;	

	//-- no cut
	int nHFplusRechitNoCut = 0;						  	
	int nHFminusRechitNoCut = 0;						  	

	for(int irh = 0; irh < nHFRechit; ++irh) {

	  if(HFRechitEta->at(irh) > 0) {
	    hHFplusRechitEnergyNoCut->Fill(HFscaling*HFRechitEnergy->at(irh),weight);
	    nHFplusRechitNoCut++;
	  }
	    
	  if(HFRechitEta->at(irh) < 0) {
	    hHFminusRechitEnergyNoCut->Fill(HFscaling*HFRechitEnergy->at(irh),weight);	
	    nHFminusRechitNoCut++;
	  }

	}
	
	//cout<<"nHFplusRechitNoCut = "<<nHFplusRechitNoCut<<endl;	
	//cout<<"nHFminusRechitNoCut = "<<nHFminusRechitNoCut<<endl;

	hnHFplusRechitNoCut->Fill(nHFplusRechitNoCut,weight);										   
	hnHFminusRechitNoCut->Fill(nHFminusRechitNoCut,weight);										   

	//-- different cuts
	for(int ihf = 0; ihf < nHFcut; ++ihf) {
	  
	  for(int irh = 0; irh < nHFRechit; ++irh) {

	    if(HFscaling*HFRechitEnergy->at(irh) > EHFcut[ihf]) {

	      //-- HF plus alone
	      if(HFRechitEta->at(irh) > 0) {
		hHFplusRechitEnergy[ihf]->Fill(HFscaling*HFRechitEnergy->at(irh),weight);		
		hHFplusRechitEt[ihf]->Fill(HFRechitEt->at(irh),weight);		
		hHFplusRechitTime[ihf]->Fill(HFRechitTime->at(irh),weight);		  
		hHFplusRechitEta[ihf]->Fill(HFRechitEta->at(irh),weight);		  
		hHFplusRechitPhi[ihf]->Fill(HFRechitPhi->at(irh),weight);		             
		hHFplusRechitEtaPhi[ihf]->Fill(HFRechitEta->at(irh),HFRechitPhi->at(irh),weight); //-- HFscaling*HFRechitEnergy->at(irh) - E weighted
	      
		hHFplusRechitDepth[ihf]->Fill(HFRechitDepth->at(irh),weight);	      
		
		if(HFRechitDepth->at(irh) == 1) {
		  hHFplusRechitiEtaDepth1[ihf]->Fill(HFRechitiEta->at(irh),weight);    
		  hHFplusRechitiPhiDepth1[ihf]->Fill(HFRechitiPhi->at(irh),weight);    
		  hHFplusRechitiEtaiPhiDepth1[ihf]->Fill(HFRechitiEta->at(irh),HFRechitiPhi->at(irh),weight); 
		}

		if(HFRechitDepth->at(irh) == 2) {
		  hHFplusRechitiEtaDepth2[ihf]->Fill(HFRechitiEta->at(irh),weight);   
		  hHFplusRechitiPhiDepth2[ihf]->Fill(HFRechitiPhi->at(irh),weight);    		  
		  hHFplusRechitiEtaiPhiDepth2[ihf]->Fill(HFRechitiEta->at(irh),HFRechitiPhi->at(irh),weight); 
		}

		nHFplusRechit_HFcut[ihf]++;					       
		HFplusRechitEtot_HFcut[ihf]+=HFscaling*HFRechitEnergy->at(irh);		       
		                                                                       
		if(HFscaling*HFRechitEnergy->at(irh) >  HFplusRechitLeadingEnergy_HFcut[ihf])   
		  HFplusRechitLeadingEnergy_HFcut[ihf] = HFscaling*HFRechitEnergy->at(irh);     
	      } //-- end HF plus alone

	      //-- HF minus alone
	      if(HFRechitEta->at(irh) < 0) {
		hHFminusRechitEnergy[ihf]->Fill(HFscaling*HFRechitEnergy->at(irh),weight);		
		hHFminusRechitEt[ihf]->Fill(HFRechitEt->at(irh),weight);		
		hHFminusRechitTime[ihf]->Fill(HFRechitTime->at(irh),weight);		  
		hHFminusRechitEta[ihf]->Fill(HFRechitEta->at(irh),weight);		  
		hHFminusRechitPhi[ihf]->Fill(HFRechitPhi->at(irh),weight);		             
		hHFminusRechitEtaPhi[ihf]->Fill(HFRechitEta->at(irh),HFRechitPhi->at(irh),weight);
	      
		hHFminusRechitDepth[ihf]->Fill(HFRechitDepth->at(irh),weight);	      
		
		if(HFRechitDepth->at(irh) == 1) {
		  hHFminusRechitiEtaDepth1[ihf]->Fill(HFRechitiEta->at(irh),weight);    
		  hHFminusRechitiPhiDepth1[ihf]->Fill(HFRechitiPhi->at(irh),weight);    
		  hHFminusRechitiEtaiPhiDepth1[ihf]->Fill(HFRechitiEta->at(irh),HFRechitiPhi->at(irh),weight); 
		}

		if(HFRechitDepth->at(irh) == 2) {
		  hHFminusRechitiEtaDepth2[ihf]->Fill(HFRechitiEta->at(irh),weight);   
		  hHFminusRechitiPhiDepth2[ihf]->Fill(HFRechitiPhi->at(irh),weight);    		  
		  hHFminusRechitiEtaiPhiDepth2[ihf]->Fill(HFRechitiEta->at(irh),HFRechitiPhi->at(irh),weight); 
		}

		nHFminusRechit_HFcut[ihf]++;					       
		HFminusRechitEtot_HFcut[ihf]+=HFscaling*HFRechitEnergy->at(irh);		       
		                                                                       
		if(HFscaling*HFRechitEnergy->at(irh) >  HFminusRechitLeadingEnergy_HFcut[ihf])   
		  HFminusRechitLeadingEnergy_HFcut[ihf] = HFscaling*HFRechitEnergy->at(irh);     
	      } //-- end HF minus alone

	    } //-- end EHFcut
	  } //-- end loop over HFRechit

	  //-- HF plus alone
	  hnHFplusRechit[ihf]->Fill(nHFplusRechit_HFcut[ihf],weight);			   
	  hHFplusRechitEtot[ihf]->Fill(HFplusRechitEtot_HFcut[ihf],weight);		   
	  if(HFplusRechitLeadingEnergy_HFcut[ihf] > 0) hHFplusRechitLeadingEnergy[ihf]->Fill(HFplusRechitLeadingEnergy_HFcut[ihf],weight);

	  //-- HF minus alone
	  hnHFminusRechit[ihf]->Fill(nHFminusRechit_HFcut[ihf],weight);			   
	  hHFminusRechitEtot[ihf]->Fill(HFminusRechitEtot_HFcut[ihf],weight);		   
	  if(HFminusRechitLeadingEnergy_HFcut[ihf] > 0) hHFminusRechitLeadingEnergy[ihf]->Fill(HFminusRechitLeadingEnergy_HFcut[ihf],weight);

	} //-- end loop over nHFcut

	//-- vertex and tracks

	double nFwd[n_delta_eta];
	for(int i = 0; i < n_delta_eta; ++i) nFwd[i] = 0;
	double nBwd[n_delta_eta];
	for(int i = 0; i < n_delta_eta; ++i) nBwd[i] = 0;
       
	double leading_track_pT = 0;
	double leading_track_eta = 0;
	double leading_track_phi = 0;

	//-- exactly one good vertex
	hBSx->Fill(BSx,weight);
	hBSy->Fill(BSy,weight);
	hBSz->Fill(BSz,weight);

	hTrackNum->Fill(Vntrack_gv[0],weight);
	
	//-- matching - monte carlo only
	if(!isData) {

	  //-- loop over generated particles		       
	  int ngenparticle = ParticleEnergy->size();
	  
	  for(int igen = 0; igen < ngenparticle; ++igen) {       
	    
	    if(ParticleStatus->at(igen) != 1) continue;	       
	    if(ParticleCharge->at(igen) == 0) continue;	       
	    if(TMath::Abs(ParticleEta->at(igen)) > 2.4) continue;
	    if(ParticlePt->at(igen) < 0.5) continue;

	    nch_matching++;

	    ChPt.push_back(ParticlePt->at(igen));
	    ChEta.push_back(ParticleEta->at(igen));
	    ChPhi.push_back(ParticlePhi->at(igen));
	    
	  }
	  
	  //-- matching - loop over tracks associated to the primary vertex - exactly one vertex	  
	  for (int itrack = 0; itrack < VnTrack->at(0); itrack++) {
	    
	    if(VTrackPt->at(0).at(itrack) < 0.5) continue; //-- cut on track pT			     
	    if(VTrackQuality->at(0).at(itrack) != 2) continue; //-- cut on track quality		     					    

	    int npixelhit =  VTrackValidPixelHit->at(0).at(itrack);
	    double abseta = TMath::Abs(VTrackEta->at(0).at(itrack));

	    if(npixelhit < 2) continue;  //-- cut on pixel hit	           
	    if(npixelhit < 3 && pixel_cut) continue;  //-- cut on pixel hit

	    if(npixelhit < 3 && abseta < 1 && combined_cut) continue;
	    if(npixelhit < 2 && abseta > 1 && combined_cut) continue;
	    ntrack_matching++;
	    
	    double dRmin = 10000;
	    int igenmin = - 1;
	    int itrackmin = -1;
	    
	    if(debugMCh) cout<<"size of the charged particle collection = "<<ChPt.size()<<endl;
	    if(debugMCh) cout<<"size of the matched charged particle collection = "<<MChPt.size()<<endl<<endl;

	    //-- try to match a charged particle ChPt, ChEta, ChPhi
	    for(unsigned int igen = 0; igen < ChPt.size(); igen++) {	    
	      double deta = VTrackEta->at(0).at(itrack) - ChEta.at(igen);
	      double dphi = VTrackPhi->at(0).at(itrack) - ChPhi.at(igen);
	      double dR2 = deta*deta + dphi*dphi;
	      double dR = TMath::Sqrt(dR2);
	      
	      if(dR > dRmin) continue;
	      dRmin = dR;
	      igenmin = igen;
	      itrackmin = itrack;
	    }

	    //-- if matching - save charged particle - remove charged particle before going to next track
	    if(dRmin > 0.05) continue;
	    
	    MChPt.push_back(ChPt.at(igenmin));
	    MChEta.push_back(ChEta.at(igenmin));
	    MChPhi.push_back(ChPhi.at(igenmin));
	    
	    ChPt.erase(ChPt.begin()+igenmin);    
	    ChEta.erase(ChEta.begin()+igenmin);  
	    ChPhi.erase(ChPhi.begin()+igenmin);  
	    
	    MTPt.push_back(VTrackPt->at(0).at(itrackmin));
	    MTEta.push_back(VTrackEta->at(0).at(itrackmin));
	    MTPhi.push_back(VTrackPhi->at(0).at(itrackmin));
	  } //-- matching - end loop over tracks assocciated to the primary vertex - exactly one vertex

	  if(debugMCh) cout<<"size of the charged particle collection = "<<ChPt.size()<<endl;
	  if(debugMCh) cout<<"size of the matched charged particle collection = "<<MChPt.size()<<endl<<endl;
	  
	  if(debugMCh) {
	    cout<<endl;
	    cout<<"N track = "<<ntrack_matching<<endl;
	    cout<<"N ch = "<<nch_matching<<endl;
	    cout<<"N matched track = "<<MTPt.size()<<endl<<endl;
	    cout<<"N matched ch = "<<MChPt.size()<<endl;
	    getchar();
	  }
	} //-- end matching - monte carlo only

	//-- loop on matched charged particles - exactly one good vertex
	for(unsigned int ich = 0; ich < MChPt.size(); ich++) {
	  hMChPt->Fill(MChPt.at(ich),weight);
	  hMChEta->Fill(MChEta.at(ich),weight);
	  hMChPhi->Fill(MChPhi.at(ich),weight);

	  hMChMTPt->Fill(MChPt.at(ich),MTPt.at(ich),weight);
	  hMChMTEta->Fill(MChEta.at(ich),MTEta.at(ich),weight);
	  hMChMTPhi->Fill(MChPhi.at(ich),MTPhi.at(ich),weight);
	}
	
	int ngoodtrack = 0;
	int nallsel = 0;

	//-- loop on tracks - exactly one good vertex
	for (int itrack = 0; itrack < VnTrack->at(0); itrack++) {	      

	  //-- cout<<itrack+1<<") quality = "<<VTrackQuality->at(0).at(itrack)<<endl;

	  if(VTrackPt->at(0).at(itrack) < 0.5) continue; //-- cut on track pT
	  hTrackQuality->Fill(VTrackQuality->at(0).at(itrack),weight);

	  if(VTrackQuality->at(0).at(itrack) != 2) continue; //-- cut on track quality

	  int npixelhit = VTrackValidPixelHit->at(0).at(itrack);
	  double abseta = TMath::Abs(VTrackEta->at(0).at(itrack));

	  if(npixelhit < 2) continue;  //-- cut on pixel hit
	  if(npixelhit < 3 && pixel_cut) continue;  //-- cut on pixel hit

	  if(npixelhit < 3 && abseta < 1 && combined_cut) continue;
	  if(npixelhit < 2 && abseta > 1 && combined_cut) continue;
	  ngoodtrack++;

	  for(int i = 0; i < n_delta_eta; ++i) {
	    if(VTrackEta->at(0).at(itrack) > eta_plus_min[i] && VTrackEta->at(0).at(itrack) < eta_plus_max[i])
	      nFwd[i]+=1;
	    if(VTrackEta->at(0).at(itrack) > eta_minus_min[i] && VTrackEta->at(0).at(itrack) < eta_minus_max[i])
	      nBwd[i]+=1;
	  }	       
	  
	  if(PUdebug) cout<<"weight = "<<weight<<endl;

	  hTrackPt->Fill(VTrackPt->at(0).at(itrack),weight);
	  hTrackEta->Fill(VTrackEta->at(0).at(itrack),weight);
	  hTrackPhi->Fill(VTrackPhi->at(0).at(itrack),weight);

	  if(TMath::Abs(VTrackEta->at(0).at(itrack)) < 1) hTrackPhiCentral->Fill(VTrackPhi->at(0).at(itrack),weight);
	  if(TMath::Abs(VTrackEta->at(0).at(itrack)) > 1) hTrackPhiFwd->Fill(VTrackPhi->at(0).at(itrack),weight);

	  hTrackEtaPhi->Fill(VTrackEta->at(0).at(itrack),VTrackPhi->at(0).at(itrack),weight);
	  hTrackEtaPhiPixel->Fill(VTrackEta->at(0).at(itrack),VTrackPhi->at(0).at(itrack),weight*VTrackValidPixelHit->at(0).at(itrack));
	  hTrackEtaPhiStrip->Fill(VTrackEta->at(0).at(itrack),VTrackPhi->at(0).at(itrack),weight*VTrackValidStripHit->at(0).at(itrack));

	  if(VTrackPt->at(0).at(itrack) > leading_track_pT) {
	    leading_track_pT = VTrackPt->at(0).at(itrack);
	    leading_track_eta = VTrackEta->at(0).at(itrack);
	    leading_track_phi = VTrackPhi->at(0).at(itrack);
	  }
	  
	  for (int ieta = 0; ieta < neta_edge-1; ieta++) {
	    if(abseta > eta_edge[ieta] && abseta < eta_edge[ieta+1]) 
	      hTrackPtEta[ieta]->Fill(VTrackPt->at(0).at(itrack),weight);
	  }
	  
	  if(VTrackPt->at(0).at(itrack) != 0) {
	    double pTerror = 100*VTrackPtError->at(0).at(itrack)/VTrackPt->at(0).at(itrack);
	    hTrackPtError->Fill(pTerror,weight);
	    hTrackPtErrorPt->Fill(VTrackPt->at(0).at(itrack),pTerror,weight);
	    hTrackPtErrorEta->Fill(VTrackEta->at(0).at(itrack),pTerror,weight);
	    hTrackPtErrorNhit->Fill(VTrackValidHit->at(0).at(itrack),pTerror,weight);
	  }
	  
	  hTrackEtaError->Fill(VTrackEtaError->at(0).at(itrack),weight);
	  hTrackPhiError->Fill(VTrackPhiError->at(0).at(itrack),weight);
	  
	  if(VTrackdxyError->at(0).at(itrack) != 0) {
	    double dxysigmaxy = VTrackdxy->at(0).at(itrack)/VTrackdxyError->at(0).at(itrack);
	    hTrackdxy->Fill(TMath::Abs(dxysigmaxy),weight);
	    hTrackdxyEta->Fill(VTrackEta->at(0).at(itrack),TMath::Abs(dxysigmaxy),weight);
	    hTrackdxyNhit->Fill(VTrackValidHit->at(0).at(itrack),TMath::Abs(dxysigmaxy),weight);
	    
	    hTrackdxyEta2D->Fill(VTrackEta->at(0).at(itrack),dxysigmaxy,weight); 
	    hTrackdxyNum2D->Fill(Vntrack_gv[0],dxysigmaxy,weight); 
	  }
	  
	  if(VTrackdzError->at(0).at(itrack) != 0) {
	    double dzsigmaz = VTrackdz->at(0).at(itrack)/VTrackdzError->at(0).at(itrack);
	    hTrackdz->Fill(TMath::Abs(dzsigmaz),weight);
	    hTrackdzEta->Fill(VTrackEta->at(0).at(itrack),TMath::Abs(dzsigmaz),weight);
	    hTrackdzNhit->Fill(VTrackValidHit->at(0).at(itrack),TMath::Abs(dzsigmaz),weight);
	    
	    hTrackdzEta2D->Fill(VTrackEta->at(0).at(itrack),dzsigmaz,weight); 
	    hTrackdzNum2D->Fill(Vntrack_gv[0],dzsigmaz,weight); 
	  }	     
	  
	  hTrackchi2->Fill(VTrackchi2->at(0).at(itrack),weight);
	  hTrackndof->Fill(VTrackndof->at(0).at(itrack),weight);
	  hTrackchi2ndof->Fill(VTrackchi2ndof->at(0).at(itrack),weight);
	  hTrackValidHit->Fill(VTrackValidHit->at(0).at(itrack),weight);
	  hTrackValidPixelHit->Fill(VTrackValidPixelHit->at(0).at(itrack),weight);
	  hTrackValidStripHit->Fill(VTrackValidStripHit->at(0).at(itrack),weight);
	  hTrackWeight->Fill(VTrackWeight->at(0).at(itrack),weight);
	  
	  hTrackValidHitEta->Fill(VTrackEta->at(0).at(itrack),VTrackValidHit->at(0).at(itrack),weight);
	  hTrackValidPixelHitEta->Fill(VTrackEta->at(0).at(itrack),VTrackValidPixelHit->at(0).at(itrack),weight);
	  hTrackValidStripHitEta->Fill(VTrackEta->at(0).at(itrack),VTrackValidStripHit->at(0).at(itrack),weight);
	  
	  //-- matching vertex track - all track
	  double DeltaRmin = 100;
	  double DeltapTmin = 100;

	  int iallmin = 0;
	  int ivertexmin = 0;

	  for (int iall = 0; iall < nTrack; iall++) {
	    
	    if(TrackPt->at(iall) < 0.5) continue; 
	    if(TrackQuality->at(iall) != 2) continue; 
	    if(TrackValidHit->at(iall) < 3) continue;  
	    	    
	    double DeltaEta = TrackEta->at(iall) - VTrackEta->at(0).at(itrack);
	    double DeltaPhi = TrackPhi->at(iall) - VTrackPhi->at(0).at(itrack);
	    double DeltaR = TMath::Sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);
	    
	    if(DeltaR > DeltaRmin) continue;

	    DeltaRmin = DeltaR;
	    iallmin = iall;
	    ivertexmin = itrack;

	  } //-- end matching vertex track - all track

	  //-- Delta R condition
	  if(DeltaRmin < 0.05) {
	    
	    nallsel++;

	    DeltapTmin = 100*(TrackPt->at(iallmin) - VTrackPt->at(0).at(ivertexmin));
	    if(VTrackPt->at(0).at(ivertexmin) > 0) DeltapTmin/=VTrackPt->at(0).at(ivertexmin);	    

	    hCheckDeltaR->Fill(DeltaRmin,weight);
	    hCheckDeltapT->Fill(DeltapTmin,weight);

	    if(HitDebug) {
	      getchar();
	      cout<<"track "<<nallsel<<" has "<<TrackHitR->at(iallmin).size()<<" hits"<<endl;
	      getchar();
	    }

	    //-- loop on track hits
	    for (unsigned int ihit = 0; ihit < TrackHitR->at(iallmin).size(); ihit++) {
	  
	      //-- find MaxLayer 
	      for(int idet = 0; idet < detsize; idet++) {
		if(TrackHitDet->at(iallmin).at(ihit) == idet+1) {
		  if(TrackHitSubdet->at(iallmin).at(ihit) > MaxLayer[idet]) {
		    MaxLayer[idet] = TrackHitSubdet->at(iallmin).at(ihit);
		  }
		}
	      }
	      
	      //-- all tracker													    
	      hHitEtaPhiAll->Fill(TrackHitEta->at(iallmin).at(ihit),TrackHitPhi->at(iallmin).at(ihit),weight);
	      hHitZRAll->Fill(TrackHitZ->at(iallmin).at(ihit),TrackHitR->at(iallmin).at(ihit),weight);
	      hHitXYAll->Fill(TrackHitX->at(iallmin).at(ihit),TrackHitY->at(iallmin).at(ihit),weight);     

	      int detbin = TrackHitDet->at(iallmin).at(ihit) - 1;
	      int layerbin = TrackHitSubdet->at(iallmin).at(ihit) - 1;

	      //-- each detector
	      hHitEtaPhiDet[detbin]->Fill(TrackHitEta->at(iallmin).at(ihit),TrackHitPhi->at(iallmin).at(ihit),weight); 
	      hHitZRDet[detbin]->Fill(TrackHitZ->at(iallmin).at(ihit),TrackHitR->at(iallmin).at(ihit),weight);     
	      hHitXYDet[detbin]->Fill(TrackHitX->at(iallmin).at(ihit),TrackHitY->at(iallmin).at(ihit),weight);     

	      //-- each layer

	      hHitEtaPhiLayer[detbin][layerbin]->Fill(TrackHitEta->at(iallmin).at(ihit),TrackHitPhi->at(iallmin).at(ihit),weight); 	
	      hHitZRLayer[detbin][layerbin]->Fill(TrackHitZ->at(iallmin).at(ihit),TrackHitR->at(iallmin).at(ihit),weight);     
	      hHitXYLayer[detbin][layerbin]->Fill(TrackHitX->at(iallmin).at(ihit),TrackHitY->at(iallmin).at(ihit),weight);          
	      
	      //-- debug section
	      if(!HitDebug) continue; 
	      cout<<"hit "<<ihit+1<<":"<<endl;
	      
	      cout<<"R = "<<TrackHitR->at(iallmin).at(ihit)<<endl;
	      cout<<"Theta = "<<TrackHitTheta->at(iallmin).at(ihit)<<endl;
	      cout<<"Phi = "<<TrackHitPhi->at(iallmin).at(ihit)<<endl;
	      
	      cout<<"Eta = "<<TrackHitEta->at(iallmin).at(ihit)<<endl;
	      
	      cout<<"X = "<<TrackHitX->at(iallmin).at(ihit)<<endl;
	      cout<<"Y = "<<TrackHitY->at(iallmin).at(ihit)<<endl;
	      cout<<"Z = "<<TrackHitZ->at(iallmin).at(ihit)<<endl;
	      
	      //-- cout<<"Det = "<<TrackHitDet->at(iallmin).at(ihit)<<endl;
	      //-- cout<<"Subdet = "<<TrackHitSubdet->at(iallmin).at(ihit)<<endl<<endl;
	      
	      //-- PXB = 1 - PXF = 2 - TIB = 3 - TID = 4 - TOB = 5 - TEC = 6
	      for(int idet = 0; idet < detsize; idet++) {
		if(TrackHitDet->at(iallmin).at(ihit) == idet+1) {
		  cout<<"subdetector = "<<DetName[idet]<<" = "<<TrackHitDet->at(iallmin).at(ihit)<<endl;
		  cout<<SubdetName[idet]<<" = "<<TrackHitSubdet->at(iallmin).at(ihit)<<endl<<endl;
		}
	      }
	      
	    } //-- end loop on track hits
	    
	    if(!MatchingDebug) continue;
	    cout<<nallsel<<") pT all track = "<<TrackPt->at(iallmin)<<" - pT vertex track = "<<VTrackPt->at(0).at(ivertexmin);
	    cout<<" - Delta R min = "<<DeltaRmin<<" - Delta pT min / pT = "<<DeltapTmin<<" %"<<endl;

	  } //-- end Delta R condition

	} //-- end loop on tracks - exactly one good vertex

	if(MatchingDebug){ 
	  cout<<endl;	  
	  cout<<"N vertex track sel = "<<ngoodtrack<<" - N all track sel = "<<nallsel<<endl;
	  cout<<endl;
	  getchar();
	}

	hCheckDeltaNum->Fill(ngoodtrack-nallsel,weight);

	//-- exactly one good vertex
	hVx->Fill(Vx->at(0),weight);
	hVxBS->Fill(VxBS[0],weight);
	hVy->Fill(Vy->at(0),weight);
	hVyBS->Fill(VyBS[0],weight);
	hVxy->Fill(Vx->at(0),Vy->at(0),weight);
	hVz->Fill(Vz->at(0),weight);
	hVrho->Fill(Vrho->at(0),weight);
	hVrhoBS->Fill(rhoBS[0],weight);
	hVchi2->Fill(Vchi2->at(0),weight);
	hVndof->Fill(Vndof->at(0),weight);
	hVchi2ndof->Fill(Vchi2ndof->at(0),weight);
	hVzError->Fill(VzError->at(0),weight);
	hVzErrorVz->Fill(Vz->at(0),VzError->at(0),weight);
	hVzErrorTrackNum->Fill(Vntrack_gv[0],VzError->at(0),weight);
	
	//-- leading track
	if(leading_track_pT > 0) {
	  hLeadingTrackPt->Fill(leading_track_pT,weight);
	  hLeadingTrackEta->Fill(leading_track_eta,weight);
	  hLeadingTrackPhi->Fill(leading_track_phi,weight);
	}
	
	//-- correlation coefficient
	for(int i = 0; i < n_delta_eta; ++i) {
	  hnFwd[i]->Fill(nFwd[i],weight);
	  hnBwd[i]->Fill(nBwd[i],weight);
	  hnFwdBwd[i]->Fill(nFwd[i]*nBwd[i],weight);
	}

      } //-- end loop on data - moca reco exactly one good vertex

      ChPt.clear();
      ChEta.clear();	     
      ChPhi.clear();	     
      
      MChPt.clear();	     
      MChEta.clear();	     
      MChPhi.clear();	     
      
      MTPt.clear(); 	     
      MTEta.clear();	     
      MTPhi.clear();       

      //-- monte carlo only
      if(!isData && moca_reco_sel) {

	int ngenparticle = ParticleEnergy->size();
	
	//-------------------
	//-- Leading particle
	//-------------------

	double PartPlusLeadingEnergy = 0;
	double PartPlusLeadingEta = 0;
	double PartPlusLeadingPhi = 0;
	int PartPlusLeadingIndex = -1;

	double PartMinusLeadingEnergy = 0;
	double PartMinusLeadingEta = 0;
	double PartMinusLeadingPhi = 0;
	int PartMinusLeadingIndex = -1;

	//-- loop over generated particles
	for(int igen = 0; igen < ngenparticle; ++igen) {
	  
	  if(ParticleStatus->at(igen) != 1) continue;
	  
	  //-- HF plus 2.853 < eta < 5.191      
	  if(ParticleEta->at(igen) > 2.853 && ParticleEta->at(igen) < 5.191) {
	    if(ParticleEnergy->at(igen) > PartPlusLeadingEnergy) {
	      PartPlusLeadingEnergy = ParticleEnergy->at(igen);
	      PartPlusLeadingEta = ParticleEta->at(igen);
	      PartPlusLeadingPhi = ParticlePhi->at(igen);
	      PartPlusLeadingIndex = igen;
	    }
	  }
	  
	  //-- HF minus -5.191 < eta < -2.853
	  if(ParticleEta->at(igen) > -5.191 && ParticleEta->at(igen) < -2.853) {
	    if(ParticleEnergy->at(igen) > PartMinusLeadingEnergy) {
	      PartMinusLeadingEnergy = ParticleEnergy->at(igen);
	      PartMinusLeadingEta = ParticleEta->at(igen);
	      PartMinusLeadingPhi = ParticlePhi->at(igen);
	      PartMinusLeadingIndex = igen;
	    }
	  }
	  
	} //-- end loop over generated particles

	//------------------------------
	//-- Closest to leading particle
	//------------------------------

	double dRmin_PartPlusClosest = 0.2011; //-- Delta eta = 0.10 - Delta phi = 10 degrees =  0.1745
	double PartPlusClosestEnergy = 0;
	double PartPlusClosestEta = 0;
	double PartPlusClosestPhi = 0;

	double dRmin_PartMinusClosest = 0.2011; //-- Delta eta = 0.10 - Delta phi = 10 degrees =  0.1745
	double PartMinusClosestEnergy = 0;
	double PartMinusClosestEta = 0;
	double PartMinusClosestPhi = 0;

	//-- loop over generated particles
	for(int igen = 0; igen < ngenparticle; ++igen) {
	  
	  if(ParticleStatus->at(igen) != 1) continue;
	  
	  //-- HF plus 2.853 < eta < 5.191      
	  if(ParticleEta->at(igen) > 2.853 && ParticleEta->at(igen) < 5.191) {
	    if(igen != PartPlusLeadingIndex) {
	
	      double dEta = ParticleEta->at(igen) - PartPlusLeadingEta;  
	      double dPhi = ParticlePhi->at(igen) - PartPlusLeadingPhi;  
	      double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);                  

	      if(dR > dRmin_PartPlusClosest) continue;

	      dRmin_PartPlusClosest = dR;
	      PartPlusClosestEnergy = ParticleEnergy->at(igen);
	      PartPlusClosestEta = ParticleEta->at(igen);
	      PartPlusClosestPhi = ParticlePhi->at(igen);
	    }
	  }

	  //-- HF minus -5.191 < eta < -2.853
	  if(ParticleEta->at(igen) > -5.191 && ParticleEta->at(igen) < -2.853) {
	    if(igen != PartMinusLeadingIndex) {
	
	      double dEta = ParticleEta->at(igen) - PartMinusLeadingEta;  
	      double dPhi = ParticlePhi->at(igen) - PartMinusLeadingPhi;  
	      double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);                  
	      
	      if(dR > dRmin_PartMinusClosest) continue;
	      
	      dRmin_PartMinusClosest = dR;
	      PartMinusClosestEnergy = ParticleEnergy->at(igen);
	      PartMinusClosestEta = ParticleEta->at(igen);
	      PartMinusClosestPhi = ParticlePhi->at(igen);
	    }
	  }

	} //-- end loop over generated particles

	if(HFleading) {

	  if(PartPlusLeadingEnergy < 5 && PartPlusClosestEnergy > 0) {
	    cout<<endl;
	    cout<<"leading particle in HF plus"<<endl;
	    cout<<"energy = "<<PartPlusLeadingEnergy<<" GeV"<<endl;
	    cout<<"eta = "<<PartPlusLeadingEta<<endl;
	    cout<<"phi = "<<PartPlusLeadingPhi<<endl;
	    cout<<endl;	  
	    
	    cout<<endl;
	    cout<<"particle closest to the leading particle in HF plus"<<endl;
	    cout<<"energy = "<<PartPlusClosestEnergy<<" GeV"<<endl;
	    cout<<"eta = "<<PartPlusClosestEta<<endl;
	    cout<<"phi = "<<PartPlusClosestPhi<<endl;
	    cout<<"dR = "<<dRmin_PartPlusClosest<<endl;
	    cout<<endl;	  
	    
	    cout<<endl;
	    cout<<"leading tower in HF plus"<<endl;
	    cout<<"energy = "<<HFplusTowerLeadingEnergy<<" GeV"<<endl;
	    cout<<"eta = "<<HFplusTowerLeadingEta<<endl;
	    cout<<"phi = "<<HFplusTowerLeadingPhi<<endl;
	    cout<<endl;	  

	    cout<<endl;
	    cout<<"Response 1 = "<<HFplusTowerLeadingEnergy/PartPlusLeadingEnergy<<endl;
	    cout<<"Response 2 = "<<HFplusTowerLeadingEnergy/(PartPlusLeadingEnergy+PartPlusClosestEnergy)<<endl;
	    cout<<endl;

	    getchar();
	  }

	  if(PartMinusLeadingEnergy < 5 && PartMinusClosestEnergy > 0) {
	    cout<<endl;
	    cout<<"leading particle in HF minus"<<endl;
	    cout<<"energy = "<<PartMinusLeadingEnergy<<" GeV"<<endl;
	    cout<<"eta = "<<PartMinusLeadingEta<<endl;
	    cout<<"phi = "<<PartMinusLeadingPhi<<endl;
	    cout<<endl;	 
	    
	    cout<<endl;
	    cout<<"particle closest to the leading particle in HF minus"<<endl;
	    cout<<"energy = "<<PartMinusClosestEnergy<<" GeV"<<endl;
	    cout<<"eta = "<<PartMinusClosestEta<<endl;
	    cout<<"phi = "<<PartMinusClosestPhi<<endl;
	    cout<<"dR = "<<dRmin_PartMinusClosest<<endl;
	    cout<<endl;	 
	    
	    cout<<endl;
	    cout<<"leading tower in HF minus"<<endl;
	    cout<<"energy = "<<HFminusTowerLeadingEnergy<<" GeV"<<endl;
	    cout<<"eta = "<<HFminusTowerLeadingEta<<endl;
	    cout<<"phi = "<<HFminusTowerLeadingPhi<<endl;
	    cout<<endl;	  
	    
	    cout<<endl;
	    cout<<"Response 1 = "<<HFminusTowerLeadingEnergy/PartMinusLeadingEnergy<<endl;
	    cout<<"Response 2 = "<<HFminusTowerLeadingEnergy/(PartMinusLeadingEnergy+PartMinusClosestEnergy)<<endl;
	    cout<<endl;

	    getchar();
	  }
	}

	//-----------------
	//-- HF correlation
	//-----------------

	//-- loop over HF Occupancy Cut
	for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {	  

	  //-- HF plus correlation 
	  if(HFplusTowerLeadingEnergy > 0 && PartPlusLeadingEnergy > EHFOccupancyCut[ihf]) {
	    hHFplus_correlation_energy[ihf]->Fill(PartPlusLeadingEnergy,HFplusTowerLeadingEnergy,weight);
	    hHFplus_correlation_eta[ihf]->Fill(PartPlusLeadingEta,HFplusTowerLeadingEta,weight);
	    hHFplus_correlation_phi[ihf]->Fill(PartPlusLeadingPhi,HFplusTowerLeadingPhi,weight);
	    Nmoca_reco_sel_HFplus_correlation[ihf]+=weight;
	  }
	
	  //-- HF minus correlation 
	  if(HFminusTowerLeadingEnergy > 0 && PartMinusLeadingEnergy > EHFOccupancyCut[ihf]) {
	    hHFminus_correlation_energy[ihf]->Fill(PartMinusLeadingEnergy,HFminusTowerLeadingEnergy,weight);
	    hHFminus_correlation_eta[ihf]->Fill(PartMinusLeadingEta,HFminusTowerLeadingEta,weight);
	    hHFminus_correlation_phi[ihf]->Fill(PartMinusLeadingPhi,HFminusTowerLeadingPhi,weight);
	    Nmoca_reco_sel_HFminus_correlation[ihf]+=weight;
	  }
	}

	//--------------
	//-- HF response
	//--------------

	if(PartPlusLeadingEnergy > 0 && HFplusTowerLeadingEnergy > 0)
	  hHFplusResponseEgen->Fill(PartPlusLeadingEnergy,HFplusTowerLeadingEnergy/PartPlusLeadingEnergy,weight);

	if(PartMinusLeadingEnergy > 0 && HFminusTowerLeadingEnergy > 0)
	  hHFminusResponseEgen->Fill(PartMinusLeadingEnergy,HFminusTowerLeadingEnergy/PartMinusLeadingEnergy,weight);

	if(PartPlusLeadingEnergy > 0 && HFplusTowerLeadingEnergy > 0) {
	  double PartPlusEnergy = PartPlusLeadingEnergy + PartPlusClosestEnergy;
	  hHFplusResponseBiasEgen->Fill(PartPlusLeadingEnergy,HFplusTowerLeadingEnergy/PartPlusEnergy,weight);
	}

	if(PartMinusLeadingEnergy > 0 && HFminusTowerLeadingEnergy > 0) {
	  double PartMinusEnergy = PartMinusLeadingEnergy + PartMinusClosestEnergy;
	  hHFminusResponseBiasEgen->Fill(PartMinusLeadingEnergy,HFminusTowerLeadingEnergy/PartMinusEnergy,weight);
	}

	for(int igen = 0; igen < nbin_response; ++igen) {
	  double Emin = igen*100.0/nbin_response;
	  double Emax = (igen+1)*100.0/nbin_response;
	  if(PartPlusLeadingEnergy > Emin && PartPlusLeadingEnergy < Emax) {
	    hHFplusResponse[igen]->Fill(HFplusTowerLeadingEnergy/PartPlusLeadingEnergy,weight);
	    hHFminusResponse[igen]->Fill(HFminusTowerLeadingEnergy/PartMinusLeadingEnergy,weight);

	    hHFResponse[igen]->Fill(HFplusTowerLeadingEnergy/PartPlusLeadingEnergy,weight);	  
	    hHFResponse[igen]->Fill(HFminusTowerLeadingEnergy/PartMinusLeadingEnergy,weight);
	  }
	}

	//-------------------------------
	//-- HF particle - tower matching
	//-------------------------------

	double dRmin_plus = 1000;
	double dEtamin_plus = 1000;
	double dPhimin_plus = 1000;
	int itower_plus = -1;

	double dRmin_minus = 1000;
	double dEtamin_minus = 1000;
	double dPhimin_minus = 1000;
	int itower_minus = -1;	

	//-- loop over tower
	for(int itower = 0; itower < nHFTower; ++itower) {

	  //-- HF plus matching
	  if(HFTowerEta->at(itower) > 0) {
	    
	    double dEta = HFTowerEta->at(itower) - PartPlusLeadingEta;  
	    double dPhi = HFTowerPhi->at(itower) - PartPlusLeadingPhi;  
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);                  

	    if(dR > dRmin_plus) continue;
	    dRmin_plus = dR;
	    dEtamin_plus = dEta;
	    dPhimin_plus = dPhi;
	    itower_plus = itower;
	  } 

	  //-- HF minus matching
	  if(HFTowerEta->at(itower) < 0) {
	    
	    double dEta = HFTowerEta->at(itower) - PartMinusLeadingEta;  
	    double dPhi = HFTowerPhi->at(itower) - PartMinusLeadingPhi;  
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);                  

	    if(dR > dRmin_minus) continue;
	    dRmin_minus = dR;
	    dEtamin_minus = dEta;
	    dPhimin_minus = dPhi;
	    itower_minus = itower;
	  } 
	} //-- end loop over tower      
	
	double HFplusTowerMatchedEnergy = 0;  	 
	double HFplusTowerMatchedEta = 0;	 
	double HFplusTowerMatchedPhi = 0;      

	if(itower_plus != -1) {
	  HFplusTowerMatchedEnergy = HFscaling*HFTowerEnergy->at(itower_plus);
	  HFplusTowerMatchedEta = HFTowerEta->at(itower_plus);
	  HFplusTowerMatchedPhi = HFTowerPhi->at(itower_plus);
	}                                                               

	double HFminusTowerMatchedEnergy = 0;  	 			
	double HFminusTowerMatchedEta = 0;	 			
	double HFminusTowerMatchedPhi = 0;      			
                                                                
	if(itower_minus != -1) {					
	  HFminusTowerMatchedEnergy = HFscaling*HFTowerEnergy->at(itower_minus);		
	  HFminusTowerMatchedEta = HFTowerEta->at(itower_minus);	
	  HFminusTowerMatchedPhi = HFTowerPhi->at(itower_minus);	
	}                                                               
	
	if(HFmatching) {
	  cout<<endl;
	  cout<<"leading particle in HF plus"<<endl;
	  cout<<"energy = "<<PartPlusLeadingEnergy<<" GeV"<<endl;
	  cout<<"eta = "<<PartPlusLeadingEta<<endl;
	  cout<<"phi = "<<PartPlusLeadingPhi<<endl;
	  cout<<endl;	  

	  cout<<endl;
	  cout<<"leading tower in HF plus"<<endl;
	  cout<<"energy = "<<HFplusTowerLeadingEnergy<<" GeV"<<endl;
	  cout<<"eta = "<<HFplusTowerLeadingEta<<endl;
	  cout<<"phi = "<<HFplusTowerLeadingPhi<<endl;
	  cout<<endl;	  

	  cout<<endl;
	  cout<<"tower matched to the leading particle in HF plus"<<endl;
	  cout<<"energy = "<<HFplusTowerMatchedEnergy<<" GeV"<<endl;
	  cout<<"eta = "<<HFplusTowerMatchedEta<<endl;
	  cout<<"phi = "<<HFplusTowerMatchedPhi<<endl;
	  cout<<endl;	

	  getchar();

	  cout<<endl;
	  cout<<"leading particle in HF minus"<<endl;
	  cout<<"energy = "<<PartMinusLeadingEnergy<<" GeV"<<endl;
	  cout<<"eta = "<<PartMinusLeadingEta<<endl;
	  cout<<"phi = "<<PartMinusLeadingPhi<<endl;
	  cout<<endl;	  

	  cout<<endl;
	  cout<<"leading tower in HF minus"<<endl;
	  cout<<"energy = "<<HFminusTowerLeadingEnergy<<" GeV"<<endl;
	  cout<<"eta = "<<HFminusTowerLeadingEta<<endl;
	  cout<<"phi = "<<HFminusTowerLeadingPhi<<endl;
	  cout<<endl;	  

	  cout<<endl;
	  cout<<"tower matched to the leading particle in HF minus"<<endl;
	  cout<<"energy = "<<HFminusTowerMatchedEnergy<<" GeV"<<endl;
	  cout<<"eta = "<<HFminusTowerMatchedEta<<endl;
	  cout<<"phi = "<<HFminusTowerMatchedPhi<<endl;
	  cout<<endl;	
	}


	//-- HF particle - tower matching													    
	hHFDeltaE->Fill(100*(HFplusTowerLeadingEnergy-HFplusTowerMatchedEnergy)/HFplusTowerLeadingEnergy,weight);	   
	hHFDeltaEta->Fill(HFplusTowerLeadingEta-HFplusTowerMatchedEta,weight);					   
	hHFDeltaPhi->Fill(HFplusTowerLeadingPhi-HFplusTowerMatchedPhi,weight);                                         

	hHFDeltaE->Fill(100*(HFminusTowerLeadingEnergy-HFminusTowerMatchedEnergy)/HFminusTowerLeadingEnergy,weight);     
	hHFDeltaEta->Fill(HFminusTowerLeadingEta-HFminusTowerMatchedEta,weight);					      
	hHFDeltaPhi->Fill(HFminusTowerLeadingPhi-HFminusTowerMatchedPhi,weight);                                         

	//-- HF plus matching
	hHFplusDeltaE->Fill(100*(HFplusTowerLeadingEnergy-HFplusTowerMatchedEnergy)/HFplusTowerLeadingEnergy,weight);
	hHFplusDeltaEta->Fill(HFplusTowerLeadingEta-HFplusTowerMatchedEta,weight);
	hHFplusDeltaPhi->Fill(HFplusTowerLeadingPhi-HFplusTowerMatchedPhi,weight);                                         

	//-- HF minus matching												   
	hHFminusDeltaE->Fill(100*(HFminusTowerLeadingEnergy-HFminusTowerMatchedEnergy)/HFminusTowerLeadingEnergy,weight);	   
	hHFminusDeltaEta->Fill(HFminusTowerLeadingEta-HFminusTowerMatchedEta,weight);					   
	hHFminusDeltaPhi->Fill(HFminusTowerLeadingPhi-HFminusTowerMatchedPhi,weight);                                         

	//---------------------	
	//-- HF tower occupancy
	//---------------------

	int ngenparticle_HFplus10 = 0;
	int ngenparticle_HFplus20 = 0;
	
	int ngenparticle_HFminus10 = 0;
	int ngenparticle_HFminus20 = 0;	
	
	//-- loop over HF Occupancy Cut
	for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
	  
	  //-- loop over generated particles
	  for(int igen = 0; igen < ngenparticle; ++igen) {
	    
	    if(ParticleStatus->at(igen) != 1) continue;
	    
	    //-- HF plus 10 degrees tower occupancy: 2.853 < eta < 4.716	    
	    if(ParticleEta->at(igen) > 2.853 && ParticleEta->at(igen) < 4.716 && ParticleEnergy->at(igen) > EHFOccupancyCut[0])
	      ngenparticle_HFplus10++;
	    
	    //-- HF plus 20 degrees tower occupancy: 4.716 < eta < 5.191      
	    if(ParticleEta->at(igen) > 4.716 && ParticleEta->at(igen) < 5.191 && ParticleEnergy->at(igen) > EHFOccupancyCut[0])
	      ngenparticle_HFplus20++;
	    
	    //-- HF minus 20 degrees tower occupancy: -5.191 < eta < -4.716 
	    if(ParticleEta->at(igen) > -5.191 && ParticleEta->at(igen) < -4.716 && ParticleEnergy->at(igen) > EHFOccupancyCut[0])
	      ngenparticle_HFminus20++;

	    //-- HF minus 10 degrees tower occupancy: -4.716 < eta < -2.853 
	    if(ParticleEta->at(igen) > -4.716 && ParticleEta->at(igen) < -2.853 && ParticleEnergy->at(igen) > EHFOccupancyCut[0])
	      ngenparticle_HFminus10++;
	    
	    //-- E HF Occupancy cut applied here
	    if(ParticleEnergy->at(igen) < EHFOccupancyCut[ihf]) continue;
	    
	    //-- HF plus 10 degrees tower occupancy: 2.853 < eta < 4.716	    
	    if(ParticleEta->at(igen) > 2.853 && ParticleEta->at(igen) < 4.716) 
	      hHFplus10Occupancy[ihf]->Fill(ParticleEta->at(igen),Rad2Deg(ParticlePhi->at(igen)),1);
	    
	    //-- HF plus 20 degrees tower occupancy: 4.716 < eta < 5.191
	    if(ParticleEta->at(igen) > 4.716 && ParticleEta->at(igen) < 5.191)
	      hHFplus20Occupancy[ihf]->Fill(ParticleEta->at(igen),Rad2Deg(ParticlePhi->at(igen)),1);
	    
	    //-- HF minus 20 degrees tower occupancy: -5.191 < eta < -4.716 
	    if(ParticleEta->at(igen) > -5.191 && ParticleEta->at(igen) < -4.716)
	      hHFminus20Occupancy[ihf]->Fill(ParticleEta->at(igen),Rad2Deg(ParticlePhi->at(igen)),1);
	    
	    //-- HF minus 10 degrees tower occupancy: -4.716 < eta < -2.853 
	    if(ParticleEta->at(igen) > -4.716 && ParticleEta->at(igen) < -2.853)
	      hHFminus10Occupancy[ihf]->Fill(ParticleEta->at(igen),Rad2Deg(ParticlePhi->at(igen)),1);
	    
	  } //-- end loop over generated particles 
	} //-- end loop over HF Occupancy Cut
	
	if(ngenparticle_HFplus10 > 0) Nmoca_reco_sel_HFplus10+=weight;
	if(ngenparticle_HFplus20 > 0) Nmoca_reco_sel_HFplus20+=weight;
	
	if(ngenparticle_HFminus10 > 0) Nmoca_reco_sel_HFminus10+=weight;
	if(ngenparticle_HFminus20 > 0) Nmoca_reco_sel_HFminus20+=weight;
      } //-- end loop monte carlo only

    } //-- end event loop
  
    delete tree;
    file->Close();
  
  } //-- end file loop

  //-- HF rechit eta - phi normalization
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    double scaling = 0;
    if(!isData) scaling = 1/Nmoca_reco_sel;
    if(isData) scaling = 1/Ndata_sel;
    
    hHFplusRechitEtaPhi[ihf]->Scale(scaling);
    hHFminusRechitEtaPhi[ihf]->Scale(scaling);

    hHFplusRechitiEtaiPhiDepth1[ihf]->Scale(scaling); 
    hHFplusRechitiEtaiPhiDepth2[ihf]->Scale(scaling); 

    hHFminusRechitiEtaiPhiDepth1[ihf]->Scale(scaling);    
    hHFminusRechitiEtaiPhiDepth2[ihf]->Scale(scaling);
  }

  //-- HF tower eta - phi normalization 
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    double scaling = 0;
    if(!isData) scaling = 1/Nmoca_reco_sel;
    if(isData) scaling = 1/Ndata_sel;

    hHFTowerEtaPhi[ihf]->Scale(scaling);
    hHFplusTowerEtaPhi[ihf]->Scale(scaling);
    hHFminusTowerEtaPhi[ihf]->Scale(scaling);
  }
  
  //-- HF tower occupancy normalization
  if(!isData) {
    for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
      hHFplus10Occupancy[ihf]->Scale(1/Nmoca_reco_sel_HFplus10);
      hHFplus20Occupancy[ihf]->Scale(1/Nmoca_reco_sel_HFplus20);
      
      hHFminus10Occupancy[ihf]->Scale(1/Nmoca_reco_sel_HFminus10);
      hHFminus20Occupancy[ihf]->Scale(1/Nmoca_reco_sel_HFminus20);    
    }
  }

  //-- HF correlation normalization
  //if(!isData) {
  //for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {	  
  //  hHFplus_correlation_energy[ihf]->Scale(100/Nmoca_reco_sel_HFplus_correlation[ihf]);
  //  hHFplus_correlation_eta[ihf]->Scale(100/Nmoca_reco_sel_HFplus_correlation[ihf]);
  //  hHFplus_correlation_phi[ihf]->Scale(100/Nmoca_reco_sel_HFplus_correlation[ihf]);
      
  //  hHFminus_correlation_energy[ihf]->Scale(100/Nmoca_reco_sel_HFminus_correlation[ihf]);
  //  hHFminus_correlation_eta[ihf]->Scale(100/Nmoca_reco_sel_HFminus_correlation[ihf]);
  //  hHFminus_correlation_phi[ihf]->Scale(100/Nmoca_reco_sel_HFminus_correlation[ihf]);
  //}
  //}

  //-- HF response gaussian fit
  TF1* fit_hfplus_response[nbin_response];
  double mean_hfplus_response[nbin_response];
  double mean_error_hfplus_response[nbin_response];

  TF1* fit_hfminus_response[nbin_response];	   
  double mean_hfminus_response[nbin_response];	   
  double mean_error_hfminus_response[nbin_response];

  TF1* fit_hf_response[nbin_response];	    
  double mean_hf_response[nbin_response];	    
  double mean_error_hf_response[nbin_response];

  if(!isData) {

    double min = 0;     
    double max = 2;     
    double DeltaM = 1.1;

    for(int igen = 0; igen < nbin_response; ++igen) {
      //-- HF plus
      fit_hfplus_response[igen] = FitGauss(hHFplusResponse[igen],min,max,DeltaM);
      
      mean_hfplus_response[igen] = fit_hfplus_response[igen]->GetParameter(1);
      mean_error_hfplus_response[igen] = fit_hfplus_response[igen]->GetParError(1);
      
      hHFplusMeanResponse->SetBinContent(igen+1,mean_hfplus_response[igen]);
      hHFplusMeanResponse->SetBinError(igen+1,mean_error_hfplus_response[igen]);     
      
      //-- HF minus
      fit_hfminus_response[igen] = FitGauss(hHFminusResponse[igen],min,max,DeltaM);		   
      
      mean_hfminus_response[igen] = fit_hfminus_response[igen]->GetParameter(1);	   
      mean_error_hfminus_response[igen] = fit_hfminus_response[igen]->GetParError(1);  
      
      hHFminusMeanResponse->SetBinContent(igen+1,mean_hfminus_response[igen]);	   
      hHFminusMeanResponse->SetBinError(igen+1,mean_error_hfminus_response[igen]);         
      
      //-- HF 									 
      fit_hf_response[igen] = FitGauss(hHFResponse[igen],min,max,DeltaM);		   	 
      
      mean_hf_response[igen] = fit_hf_response[igen]->GetParameter(1);	   	 
      mean_error_hf_response[igen] = fit_hf_response[igen]->GetParError(1);  	 
      
      hHFMeanResponse->SetBinContent(igen+1,mean_hf_response[igen]);	   	 
      hHFMeanResponse->SetBinError(igen+1,mean_error_hf_response[igen]);         
    }
  }
 
  //-- Check Histo
  cout<<endl;

  //-- generated level information (no pT cut, no reweighting)										     
  if(!isData){
    CheckHisto(hChPt);
    CheckHisto(hChEta);
    CheckHisto(hChPhi);
    
    CheckHisto(hChLeadingPt);
    CheckHisto(hChLeadingEta);
    CheckHisto(hChLeadingPhi);
    CheckHisto(hChNum);                                          

    CheckHisto(h500ChPt);						 
    CheckHisto(h500ChEta);						 
    CheckHisto(h500ChPhi);						 
    								 
    CheckHisto(h500ChLeadingPt);					 
    CheckHisto(h500ChLeadingEta);					 
    CheckHisto(h500ChLeadingPhi);					 
    CheckHisto(h500ChNum);                                       

    CheckHisto(hHadPt);						 
    CheckHisto(hHadEta);						 
    CheckHisto(hHadPhi);						 
    CheckHisto(hHadNum);                                          
  }

  //-- selection
  CheckHisto(hselection_data);
  CheckHisto(hselection_moca_reco);

  //-- beam spot
  CheckHisto(hBSx);
  CheckHisto(hBSy);
  CheckHisto(hBSz);

  //-- HF tower occupancy
  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    CheckHisto(hHFplus10Occupancy[ihf]);
    CheckHisto(hHFplus20Occupancy[ihf]);

    CheckHisto(hHFminus10Occupancy[ihf]);
    CheckHisto(hHFminus20Occupancy[ihf]);
  }

  //-- HF correlation
  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    CheckHisto(hHFplus_correlation_energy[ihf]);
    CheckHisto(hHFplus_correlation_eta[ihf]);
    CheckHisto(hHFplus_correlation_phi[ihf]);
    
    CheckHisto(hHFminus_correlation_energy[ihf]);
    CheckHisto(hHFminus_correlation_eta[ihf]);
    CheckHisto(hHFminus_correlation_phi[ihf]);
  }
  
  //-- HF response
  CheckHisto(hHFplusResponseEgen);
  CheckHisto(hHFminusResponseEgen);

  CheckHisto(hHFplusResponseBiasEgen); 
  CheckHisto(hHFminusResponseBiasEgen);

  for(int igen = 0; igen < nbin_response; ++igen){ 
    CheckHisto(hHFResponse[igen]);
    CheckHisto(hHFplusResponse[igen]);
    CheckHisto(hHFminusResponse[igen]);
  }

  CheckHisto(hHFMeanResponse);
  CheckHisto(hHFplusMeanResponse);
  CheckHisto(hHFminusMeanResponse);

  //-- HF matching
  CheckHisto(hHFDeltaE);  
  CheckHisto(hHFDeltaEta);
  CheckHisto(hHFDeltaPhi);

  CheckHisto(hHFplusDeltaE);
  CheckHisto(hHFplusDeltaEta);
  CheckHisto(hHFplusDeltaPhi);

  CheckHisto(hHFminusDeltaE);  
  CheckHisto(hHFminusDeltaEta);
  CheckHisto(hHFminusDeltaPhi);

  //-- correlation coefficient
  for(int i = 0; i < n_delta_eta; ++i) {
    CheckHisto(hnFwd[i]);
    CheckHisto(hnBwd[i]);
    CheckHisto(hnFwdBwd[i]);
  }                                         

  //-- HF rechit	
  
  CheckHisto(hHFplusRechitEnergyNoCut);	
  CheckHisto(hHFminusRechitEnergyNoCut);

  CheckHisto(hnHFplusRechitNoCut);	
  CheckHisto(hnHFminusRechitNoCut);
		    
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    CheckHisto(hnHFplusRechit[ihf]);		    
    CheckHisto(hHFplusRechitEtot[ihf]);	    
    CheckHisto(hHFplusRechitLeadingEnergy[ihf]); 
                             	    
    CheckHisto(hHFplusRechitEnergy[ihf]);	    
    CheckHisto(hHFplusRechitEt[ihf]);	    
    CheckHisto(hHFplusRechitTime[ihf]);  	    
    CheckHisto(hHFplusRechitEta[ihf]);	    
    CheckHisto(hHFplusRechitPhi[ihf]);           
    CheckHisto(hHFplusRechitEtaPhi[ihf]);	    
    				    
    CheckHisto(hHFplusRechitDepth[ihf]);	    
                                    
    CheckHisto(hHFplusRechitiEtaDepth1[ihf]);    
    CheckHisto(hHFplusRechitiEtaDepth2[ihf]);    
    
    CheckHisto(hHFplusRechitiPhiDepth1[ihf]);    
    CheckHisto(hHFplusRechitiPhiDepth2[ihf]);    
    
    CheckHisto(hHFplusRechitiEtaiPhiDepth1[ihf]);
    CheckHisto(hHFplusRechitiEtaiPhiDepth2[ihf]);

    CheckHisto(hnHFminusRechit[ihf]);		    
    CheckHisto(hHFminusRechitEtot[ihf]);	    
    CheckHisto(hHFminusRechitLeadingEnergy[ihf]); 
                             	    
    CheckHisto(hHFminusRechitEnergy[ihf]);	    
    CheckHisto(hHFminusRechitEt[ihf]);	    
    CheckHisto(hHFminusRechitTime[ihf]);  	    
    CheckHisto(hHFminusRechitEta[ihf]);	    
    CheckHisto(hHFminusRechitPhi[ihf]);           
    CheckHisto(hHFminusRechitEtaPhi[ihf]);	    
    				    
    CheckHisto(hHFminusRechitDepth[ihf]);	    
                                    
    CheckHisto(hHFminusRechitiEtaDepth1[ihf]);    
    CheckHisto(hHFminusRechitiEtaDepth2[ihf]);    
    
    CheckHisto(hHFminusRechitiPhiDepth1[ihf]);    
    CheckHisto(hHFminusRechitiPhiDepth2[ihf]);    
    
    CheckHisto(hHFminusRechitiEtaiPhiDepth1[ihf]);
    CheckHisto(hHFminusRechitiEtaiPhiDepth2[ihf]);
  }                                                      
  cout<<endl;

  //-- HF tower 
  
  CheckHisto(hHFplusTowerEnergyNoCut);
  CheckHisto(hHFminusTowerEnergyNoCut);

  CheckHisto(hnHFplusTowerNoCut); 
  CheckHisto(hnHFminusTowerNoCut);

  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    CheckHisto(hnHFTower[ihf]);
    CheckHisto(hHFTowerEtot[ihf]);
    CheckHisto(hHFTowerEnergy[ihf]);
    CheckHisto(hHFTowerEta[ihf]);
    CheckHisto(hHFTowerPhi[ihf]);
    CheckHisto(hHFTowerEtaPhi[ihf]);
    CheckHisto(hHFTowerLeadingEnergy[ihf]);
  }
  
  cout<<endl;
  
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    CheckHisto(hnHFplusTower[ihf]);		   
    CheckHisto(hHFplusTowerEtot[ihf]);	   
    CheckHisto(hHFplusTowerEnergy[ihf]);	   
    CheckHisto(hHFplusTowerEta[ihf]);	   
    CheckHisto(hHFplusTowerPhi[ihf]);    
    CheckHisto(hHFplusTowerEtaPhi[ihf]);
    CheckHisto(hHFplusTowerLeadingEnergy[ihf]);
  }
  
  cout<<endl;

  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    CheckHisto(hnHFminusTower[ihf]);	       
    CheckHisto(hHFminusTowerEtot[ihf]);	       
    CheckHisto(hHFminusTowerEnergy[ihf]);       
    CheckHisto(hHFminusTowerEta[ihf]);	       
    CheckHisto(hHFminusTowerPhi[ihf]);    
    CheckHisto(hHFminusTowerEtaPhi[ihf]);
    CheckHisto(hHFminusTowerLeadingEnergy[ihf]);
  }
  
  cout<<endl;

  //-- matched charged particle distributions - exactly one good vertex									 
  CheckHisto(hMChPt);
  CheckHisto(hMChEta);
  CheckHisto(hMChPhi);                                                   

  CheckHisto(hMChMTPt);
  CheckHisto(hMChMTEta);
  CheckHisto(hMChMTPhi);

  //-- track distributions exactly one good vertex
  CheckHisto(hTrackNum);
  CheckHisto(hTrackQuality);
  CheckHisto(hTrackPt);
  CheckHisto(hTrackEta);

  CheckHisto(hTrackPhi);
  CheckHisto(hTrackPhiCentral);
  CheckHisto(hTrackPhiFwd);

  CheckHisto(hTrackEtaPhi);
  CheckHisto(hTrackEtaPhiPixel);
  CheckHisto(hTrackEtaPhiStrip);

  CheckHisto(hLeadingTrackPt);
  CheckHisto(hLeadingTrackEta);
  CheckHisto(hLeadingTrackPhi);

  for (int ieta = 0; ieta < neta_edge-1; ieta++) CheckHisto(hTrackPtEta[ieta]);
  CheckHisto(hTrackPtError);	
  CheckHisto(hTrackEtaError);
  CheckHisto(hTrackPhiError);
  CheckHisto(hTrackPtErrorPt);	
  CheckHisto(hTrackPtErrorEta);
  CheckHisto(hTrackPtErrorNhit);
  CheckHisto(hTrackdxy);
  CheckHisto(hTrackdz);
  CheckHisto(hTrackdxyEta);
  CheckHisto(hTrackdzEta);  
  CheckHisto(hTrackdxyNhit); 
  CheckHisto(hTrackdzNhit);  
  CheckHisto(hTrackdxyEta2D);
  CheckHisto(hTrackdzEta2D);  
  CheckHisto(hTrackdxyNum2D); 
  CheckHisto(hTrackdzNum2D);  
  CheckHisto(hTrackchi2);
  CheckHisto(hTrackndof);
  CheckHisto(hTrackchi2ndof);
  CheckHisto(hTrackValidHit);
  CheckHisto(hTrackValidPixelHit);
  CheckHisto(hTrackValidStripHit);
  CheckHisto(hTrackWeight);  

  CheckHisto(hTrackValidHitEta);
  CheckHisto(hTrackValidPixelHitEta);
  CheckHisto(hTrackValidStripHitEta);

  //-- track hit
  CheckHisto(hCheckDeltaR);
  CheckHisto(hCheckDeltapT);
  CheckHisto(hCheckDeltaNum);

  //-- all tracker													    
  CheckHisto(hHitEtaPhiAll);
  CheckHisto(hHitZRAll);	    
  CheckHisto(hHitXYAll);

  //-- each detector												      
  for (int idet = 0; idet < detsize; idet++){	
    CheckHisto(hHitEtaPhiDet[idet]);
    CheckHisto(hHitZRDet[idet]);
    CheckHisto(hHitXYDet[idet]);
    
    //-- each layer													
    for (int ilayer = 0; ilayer < nlayer[idet]; ilayer++) {  
      CheckHisto(hHitEtaPhiLayer[idet][ilayer]);
      CheckHisto(hHitZRLayer[idet][ilayer]);
      CheckHisto(hHitXYLayer[idet][ilayer]);
    }						     
  }                                                          
  
  cout<<endl;
  
  //-- pileup information			    
  if(!isData) {	
    CheckHisto(hnbxint);			    
    CheckHisto(hnbxearly);
    CheckHisto(hnbxlate);  
    CheckHisto(hnbxtot);

    CheckHisto(hPUint);
    CheckHisto(hPUearly);
    CheckHisto(hPUlate);
    CheckHisto(hPUtot);

    CheckHisto(hNumInt);

    CheckHisto(hPUintVNum);
    CheckHisto(hPUearlyVNum);
    CheckHisto(hPUlateVNum);
    CheckHisto(hPUtotVNum);    

    CheckHisto(hPUintVNum0);
  }                                           
  
  //-- all vertices 
  CheckHisto(hVNum);
  CheckHisto(hVNumGood);
  CheckHisto(hVNum2D);

  //-- vertex distributions exactly one good vertex
  cout<<endl;
  CheckHisto(hVx);
  CheckHisto(hVxBS);
  CheckHisto(hVy);
  CheckHisto(hVyBS);
  CheckHisto(hVxy);
  CheckHisto(hVz);
  CheckHisto(hVrho);
  CheckHisto(hVrhoBS);
  CheckHisto(hVchi2);
  CheckHisto(hVndof);
  CheckHisto(hVchi2ndof);
  CheckHisto(hVzError);
  CheckHisto(hVzErrorVz);
  CheckHisto(hVzErrorTrackNum);

  //-- primary vertices
  cout<<endl;
  CheckHisto(hPVx);
  CheckHisto(hPVy);
  CheckHisto(hPVz);
  CheckHisto(hPVrho);
  CheckHisto(hPVchi2);
  CheckHisto(hPVndof);
  CheckHisto(hPVchi2ndof);
  CheckHisto(hPVnTrack);
  CheckHisto(hPVfTrack);
  CheckHisto(hPVTrackPt);
  CheckHisto(hPVTrackEta);
  
  //-- second vertices
  cout<<endl;
  CheckHisto(hSVx);
  CheckHisto(hSVy);
  CheckHisto(hSVz);
  CheckHisto(hSVrho);
  CheckHisto(hSVchi2);
  CheckHisto(hSVndof);
  CheckHisto(hSVchi2ndof);
  CheckHisto(hSVnTrack);
  CheckHisto(hSVfTrack);

  //-- vertices correlation
  CheckHisto(hPVzSVz);
  CheckHisto(hDeltaPVzSVz);
  CheckHisto(hSVnTrackDeltaPVzSVz);
  CheckHisto(hPVnTrackSVnTrack);
  CheckHisto(hPVfTrackSVfTrack);

  cout<<endl<<endl;

  //-- write histo to output root file 
  
  //-- correlation coefficient
  double mean_fb[n_delta_eta];
  double mean_f[n_delta_eta];
  double mean_b[n_delta_eta];
  double sigma_f[n_delta_eta];
  double sigma_b[n_delta_eta];  

  double mean_fb_error[n_delta_eta]; 
  double mean_f_error[n_delta_eta];	 
  double mean_b_error[n_delta_eta];	 
  double sigma_f_error[n_delta_eta];	 
  double sigma_b_error[n_delta_eta];  
  
  double rho[n_delta_eta];
  double rho_error[n_delta_eta];

  double delta_eta[n_delta_eta];
  double delta_eta_error[n_delta_eta];

  for(int i = 0; i < n_delta_eta; ++i) {    

    mean_fb[i] = hnFwdBwd[i]->GetMean();
    mean_f[i] = hnFwd[i]->GetMean();
    mean_b[i] = hnBwd[i]->GetMean();

    sigma_f[i] = hnFwd[i]->GetRMS();
    sigma_b[i] = hnBwd[i]->GetRMS();        

    mean_fb_error[i] = hnFwdBwd[i]->GetMeanError();    
    mean_f_error[i] = hnFwd[i]->GetMeanError();	     
    mean_b_error[i] = hnBwd[i]->GetMeanError();	     
                                         
    sigma_f_error[i] = hnFwd[i]->GetRMSError();	     
    sigma_b_error[i] = hnBwd[i]->GetRMSError();        

    rho[i] = (mean_fb[i] - mean_f[i]*mean_b[i])/(sigma_f[i]*sigma_b[i]);
    
    double a = mean_fb_error[i]/(sigma_f[i]*sigma_b[i]);
    double a2 = a*a;
    double b = (mean_f_error[i]*mean_b[i])/(sigma_f[i]*sigma_b[i]);
    double b2 = b*b;
    double c = (mean_b_error[i]*mean_f[i])/(sigma_f[i]*sigma_b[i]);
    double c2 = c*c;
    double d = sigma_f_error[i]*rho[i]/sigma_f[i];
    double d2 = d*d;
    double e = sigma_b_error[i]*rho[i]/sigma_b[i];
    double e2 = e*e;

    rho_error[i] = TMath::Sqrt(a2 + b2 + c2 + d2 + e2); 

    delta_eta[i] = eta_plus_min[i]-eta_minus_max[i];
    cout<<"delta eta = "<<delta_eta[i]<<": rho = "<<rho[i]<<" +/- "<<rho_error[i]<<endl;
    delta_eta_error[i] = 0.5*delta_eta_0;
  }                            
  
  cout<<endl<<endl;
  
  if(!isData) {
    for(int igen = 0; igen < nbin_response; ++igen)
      cout<<igen+1<<") HF plus mean response = "<<mean_hfplus_response[igen]<<" +/- "<<mean_error_hfplus_response[igen]<<endl;
  
    cout<<endl<<endl;
    
    for(int igen = 0; igen < nbin_response; ++igen) 
      cout<<igen+1<<") HF minus mean response = "<<mean_hfminus_response[igen]<<" +/- "<<mean_error_hfminus_response[igen]<<endl; 
    
    cout<<endl<<endl;
    
    for(int igen = 0; igen < nbin_response; ++igen) 
      cout<<igen+1<<") HF mean response = "<<mean_hf_response[igen]<<" +/- "<<mean_error_hf_response[igen]<<endl; 
  
    cout<<endl<<endl; 
  }

  Char_t outputfile_name[200];
  sprintf(outputfile_name,"%s",file_name.Data());
  TFile* foutput = new TFile(outputfile_name,"RECREATE");
  foutput->cd();

  //-- generated level information (no pT cut, no reweighting)    
  if(!isData) {
    hChPt->Write();					       
    hChEta->Write();					       
    hChPhi->Write();					       
    
    hChLeadingPt->Write();				       
    hChLeadingEta->Write();				       
    hChLeadingPhi->Write();				       
    hChNum->Write();                                          

    h500ChPt->Write();					      
    h500ChEta->Write();					      
    h500ChPhi->Write();					      
    							      
    h500ChLeadingPt->Write();				      
    h500ChLeadingEta->Write();				      
    h500ChLeadingPhi->Write();				      
    h500ChNum->Write();                                          
    
    hHadPt->Write();						  
    hHadEta->Write();					  
    hHadPhi->Write();					  
    hHadNum->Write();                                          
  }

  //-- selection		   
  hselection_data->Write();	   
  hselection_moca_reco->Write();

  //-- beam spot   
  hBSx->Write();
  hBSy->Write();
  hBSz->Write();

  //-- HF tower occupancy
  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    hHFplus10Occupancy[ihf]->Write();
    hHFplus20Occupancy[ihf]->Write();

    hHFminus10Occupancy[ihf]->Write();
    hHFminus20Occupancy[ihf]->Write();
  }

  //-- HF correlation
  for(int ihf = 0; ihf < nHFOccupancyCut; ++ihf) {
    hHFplus_correlation_energy[ihf]->Write();
    hHFplus_correlation_eta[ihf]->Write();
    hHFplus_correlation_phi[ihf]->Write();
    
    hHFminus_correlation_energy[ihf]->Write();
    hHFminus_correlation_eta[ihf]->Write();
    hHFminus_correlation_phi[ihf]->Write();
  }
  
  //-- HF response
  hHFplusResponseEgen->Write();
  hHFminusResponseEgen->Write();

  hHFplusResponseBiasEgen->Write(); 
  hHFminusResponseBiasEgen->Write();

  for(int igen = 0; igen < nbin_response; ++igen) {
    hHFResponse[igen]->Write();
    hHFplusResponse[igen]->Write();
    hHFminusResponse[igen]->Write();
  }

  hHFMeanResponse->Write();
  hHFplusMeanResponse->Write();
  hHFminusMeanResponse->Write();

  //-- HF matching
  hHFDeltaE->Write();  
  hHFDeltaEta->Write();
  hHFDeltaPhi->Write();

  hHFplusDeltaE->Write();
  hHFplusDeltaEta->Write();
  hHFplusDeltaPhi->Write();

  hHFminusDeltaE->Write();  
  hHFminusDeltaEta->Write();
  hHFminusDeltaPhi->Write();

  //-- correlation coefficient		    
   for(int i = 0; i < n_delta_eta; ++i) {   
    hnFwd[i]->Write();		    	    
    hnBwd[i]->Write();		    	    
    hnFwdBwd[i]->Write();		    
  }                                         

   TGraphErrors* grho = new TGraphErrors(n_delta_eta,delta_eta,rho,delta_eta_error,rho_error);     
   grho->SetName("correlation-coefficient");						       
   grho->SetTitle("correlation-coefficient");						       
   grho->Write();                                                                               
   
   //-- HF rechit
   
   hHFplusRechitEnergyNoCut->Write();	
   hHFminusRechitEnergyNoCut->Write();

   hnHFplusRechitNoCut->Write(); 
   hnHFminusRechitNoCut->Write();
   
   for(int ihf = 0; ihf < nHFcut; ++ihf) {		 
     hnHFplusRechit[ihf]->Write();		    	 
     hHFplusRechitEtot[ihf]->Write();	    		 
     hHFplusRechitLeadingEnergy[ihf]->Write(); 	 
     
     hHFplusRechitEnergy[ihf]->Write();	    	 
     hHFplusRechitEt[ihf]->Write();	    		 
     hHFplusRechitTime[ihf]->Write();  	    	 
     hHFplusRechitEta[ihf]->Write();	    		 
     hHFplusRechitPhi[ihf]->Write();           	 
     hHFplusRechitEtaPhi[ihf]->Write();	    	 
     
     hHFplusRechitDepth[ihf]->Write();	    	 
     
     hHFplusRechitiEtaDepth1[ihf]->Write();    	 
     hHFplusRechitiEtaDepth2[ihf]->Write();    	 
     
     hHFplusRechitiPhiDepth1[ihf]->Write();    	 
     hHFplusRechitiPhiDepth2[ihf]->Write();    	 
     
     hHFplusRechitiEtaiPhiDepth1[ihf]->Write();	 
     hHFplusRechitiEtaiPhiDepth2[ihf]->Write();	 

     hnHFminusRechit[ihf]->Write();		    	 
     hHFminusRechitEtot[ihf]->Write();	    		 
     hHFminusRechitLeadingEnergy[ihf]->Write(); 	 
     
     hHFminusRechitEnergy[ihf]->Write();	    	 
     hHFminusRechitEt[ihf]->Write();	    		 
     hHFminusRechitTime[ihf]->Write();  	    	 
     hHFminusRechitEta[ihf]->Write();	    		 
     hHFminusRechitPhi[ihf]->Write();           	 
     hHFminusRechitEtaPhi[ihf]->Write();	    	 
     
     hHFminusRechitDepth[ihf]->Write();	    	 
     
     hHFminusRechitiEtaDepth1[ihf]->Write();    	 
     hHFminusRechitiEtaDepth2[ihf]->Write();    	 
     
     hHFminusRechitiPhiDepth1[ihf]->Write();    	 
     hHFminusRechitiPhiDepth2[ihf]->Write();    	 
     
     hHFminusRechitiEtaiPhiDepth1[ihf]->Write();	 
     hHFminusRechitiEtaiPhiDepth2[ihf]->Write();	 
   }                                                      
   
   //-- HF tower 
   
   hHFplusTowerEnergyNoCut->Write(); 
   hHFminusTowerEnergyNoCut->Write();

   hnHFplusTowerNoCut->Write(); 
   hnHFminusTowerNoCut->Write();
   
   for(int ihf = 0; ihf < nHFcut; ++ihf) {
     hnHFTower[ihf]->Write();
     hHFTowerEtot[ihf]->Write();
     hHFTowerEnergy[ihf]->Write();
     hHFTowerEta[ihf]->Write();
     hHFTowerPhi[ihf]->Write();    
     hHFTowerEtaPhi[ihf]->Write();
     hHFTowerLeadingEnergy[ihf]->Write();
     
     hnHFplusTower[ihf]->Write();		
     hHFplusTowerEtot[ihf]->Write();	   	
     hHFplusTowerEnergy[ihf]->Write();	
     hHFplusTowerEta[ihf]->Write();	   	
     hHFplusTowerPhi[ihf]->Write();   
     hHFplusTowerEtaPhi[ihf]->Write();
     hHFplusTowerLeadingEnergy[ihf]->Write();	
     
     hnHFminusTower[ihf]->Write();	       	
     hHFminusTowerEtot[ihf]->Write();	       	
     hHFminusTowerEnergy[ihf]->Write();       
     hHFminusTowerEta[ihf]->Write();	       	
     hHFminusTowerPhi[ihf]->Write();   
     hHFminusTowerEtaPhi[ihf]->Write();
     hHFminusTowerLeadingEnergy[ihf]->Write();
   }

   //-- matched charged particle distributions - exactly one good vertex	 
   hMChPt->Write();							 
   hMChEta->Write(); 							 
   hMChPhi->Write();                                                   
   
   hMChMTPt->Write();	
   hMChMTEta->Write();
   hMChMTPhi->Write();
   
   //-- track distributions exactly one good vertex
   hTrackNum->Write(); 
   hTrackQuality->Write();
   hTrackPt->Write();
   hTrackEta->Write();
   
   hTrackPhi->Write();
   hTrackPhiCentral->Write();
   hTrackPhiFwd->Write();
   
   hTrackEtaPhi->Write();	
   hTrackEtaPhiPixel->Write();
   hTrackEtaPhiStrip->Write();
   
   hLeadingTrackPt->Write(); 
   hLeadingTrackEta->Write();
   hLeadingTrackPhi->Write();
   for (int ieta = 0; ieta < neta_edge-1; ieta++) hTrackPtEta[ieta]->Write();
   hTrackPtError->Write(); 
   hTrackEtaError->Write();
   hTrackPhiError->Write();
   hTrackPtErrorPt->Write();	  
   hTrackPtErrorEta->Write();
   hTrackPtErrorNhit->Write();	
   hTrackdxy->Write();
   hTrackdz->Write();  
   hTrackdxyEta->Write(); 
   hTrackdzEta->Write();
   hTrackdxyNhit->Write(); 
   hTrackdzNhit->Write();  
   hTrackdxyEta2D->Write(); 
   hTrackdzEta2D->Write();  
   hTrackdxyNum2D->Write(); 
   hTrackdzNum2D->Write();  
   hTrackchi2->Write();    
   hTrackndof->Write();    
   hTrackchi2ndof->Write();
   hTrackValidHit->Write();
   hTrackValidPixelHit->Write();
   hTrackValidStripHit->Write();
   hTrackWeight->Write();       
   
   hTrackValidHitEta->Write();     
   hTrackValidPixelHitEta->Write();
   hTrackValidStripHitEta->Write();
   
   //-- track hit	    
   hCheckDeltaR->Write(); 
   hCheckDeltapT->Write();
   hCheckDeltaNum->Write();
   
   //-- all tracker													    
   hHitEtaPhiAll->Write();
   hHitZRAll->Write();	    
   hHitXYAll->Write();
   
   //-- each detector												      
   for (int idet = 0; idet < detsize; idet++){	
     hHitEtaPhiDet[idet]->Write();
     hHitZRDet[idet]->Write();
     hHitXYDet[idet]->Write();
     
     //-- each layer													
     for (int ilayer = 0; ilayer < nlayer[idet]; ilayer++) {  
       hHitEtaPhiLayer[idet][ilayer]->Write();
       hHitZRLayer[idet][ilayer]->Write();
       hHitXYLayer[idet][ilayer]->Write();
     }						     
   }                                                          
   
   //-- pileup information     
   if(!isData) {		      
     hnbxint->Write();   
     hnbxearly->Write(); 
     hnbxlate->Write();  
     hnbxtot->Write();	      
     
     hPUint->Write();	      
     hPUearly->Write();     
     hPUlate->Write();      
     hPUtot->Write();	      
     
     hNumInt->Write();
     
     hPUintVNum->Write();    
     hPUearlyVNum->Write();      
     hPUlateVNum->Write();   
     hPUtotVNum->Write();    
     
     hPUintVNum0->Write();    
   }                           
   
   //-- all vertices
   hVNum->Write();
   hVNumGood->Write();
   hVNum2D->Write();
   
   //-- vertex distributions exactly one good vertex
   cout<<endl;
   hVx->Write();
   hVxBS->Write();
   hVy->Write();
   hVyBS->Write();
   hVxy->Write();
   hVz->Write();
   hVrho->Write();
   hVrhoBS->Write();
   hVchi2->Write();
   hVndof->Write();
   hVchi2ndof->Write();
   hVzError->Write();
   hVzErrorVz->Write();
   hVzErrorTrackNum->Write();
   
   //-- primary vertices
   cout<<endl;
   hPVx->Write();
   hPVy->Write();
   hPVz->Write();
   hPVrho->Write();
   hPVchi2->Write();
   hPVndof->Write();
   hPVchi2ndof->Write();
   hPVnTrack->Write();
   hPVfTrack->Write();
   hPVTrackPt->Write();   
   hPVTrackEta->Write();
   
   //-- second vertices
   cout<<endl;
   hSVx->Write();
   hSVy->Write();
   hSVz->Write();
   hSVrho->Write();
   hSVchi2->Write();
   hSVndof->Write();
   hSVchi2ndof->Write();
   hSVnTrack->Write();
   hSVfTrack->Write();
   
   //-- vertices correlation
   hPVzSVz->Write();
   hDeltaPVzSVz->Write();
   hSVnTrackDeltaPVzSVz->Write();
   hPVnTrackSVnTrack->Write();
   hPVfTrackSVfTrack->Write();
   
   foutput->Close();
   cout<<"file "<<outputfile_name<<" created"<<endl;
   
   //-- some print out
   if(isData) {
     cout<<endl<<"data selection: "<<endl;
     cout<<"number of events before selection: "<<hselection_data->GetBinContent(1)<<endl;
     cout<<"number of events after run and lumi section: "<<hselection_data->GetBinContent(2)<<endl;
     cout<<"number of events after L1TT selection: "<<hselection_data->GetBinContent(3)<<endl;
     cout<<"number of events after HLT: "<<hselection_data->GetBinContent(4)<<endl;
     cout<<"number of events after at least one good vertex: "<<hselection_data->GetBinContent(5)<<endl;
     cout<<"number of events after exactly one good vertex: "<<hselection_data->GetBinContent(6)<<endl;
   }
   
   if(!isData) {
     cout<<endl<<"moca at reco level: "<<endl;
     cout<<"number of events before selection: "<<hselection_moca_reco->GetBinContent(1)<<endl;
     cout<<"number of events after at least one good vertex: "<<hselection_moca_reco->GetBinContent(2)<<endl;
     cout<<"number of events after exactly one good vertex: "<<hselection_moca_reco->GetBinContent(3)<<endl;
   }
   
   
   cout<<endl<<endl;
   cout<<"run on "<<file_nb<<" file(s)"<<endl<<endl;
   
   cout<<"total number of events - global weight: "<<Nevt_tot<<endl;
   cout<<"total number of events - no weight: "<<Nevt_tot_no_weight<<endl;
   cout<<"total number of events - only vertex weight: "<<Nevt_tot_no_pu_weight<<endl;
   cout<<"total number of events - only pu weight: "<<Nevt_tot_no_vz_weight<<endl<<endl; 
   
   long double ratio_evt_tot_1 = Nevt_tot_no_weight/Nevt_tot;
   cout<<"no weight / global weight: "<<ratio_evt_tot_1<<endl;       
   
   long double ratio_evt_tot_2 = Nevt_tot_no_weight/Nevt_tot_no_pu_weight;
   cout<<"no weight / only vertex weight: "<<ratio_evt_tot_2<<endl;       
   
   long double ratio_evt_tot_3 = Nevt_tot_no_weight/Nevt_tot_no_vz_weight;
   cout<<"no weight / only pu weight: "<<ratio_evt_tot_3<<endl<<endl;       
   
   if(isData) cout<<"selected events in data - global weight: "<<Ndata_sel<<endl;
   if(isData) cout<<"selected events in data - no weight: "<<Ndata_sel_no_weight<<endl; 
   if(isData) cout<<"selected events in data - only vertex weight: "<<Ndata_sel_no_pu_weight<<endl;
   if(isData) cout<<"selected events in data - only pu weight: "<<Ndata_sel_no_vz_weight<<endl<<endl;
   
   long double ratio_data_sel_1 = 0;
   if(isData)  ratio_data_sel_1 = Ndata_sel_no_weight/Ndata_sel;		   
   if(isData) cout<<"no weight / global weight: "<<ratio_data_sel_1<<endl;       	   
   
   long double ratio_data_sel_2 = 0;
   if(isData)  ratio_data_sel_2 = Ndata_sel_no_weight/Ndata_sel_no_pu_weight;  
   if(isData) cout<<"no weight / only vertex weight: "<<ratio_data_sel_2<<endl;         
   
   long double ratio_data_sel_3 = 0;
   if(isData)  ratio_data_sel_3 = Ndata_sel_no_weight/Ndata_sel_no_vz_weight;  
   if(isData) cout<<"no weight / only pu weight: "<<ratio_data_sel_3<<endl<<endl;       
   
   if(!isData) cout<<"selected events in mc at reco level - global weight: "<<Nmoca_reco_sel<<endl;
   if(!isData) cout<<"selected events in mc at reco level - no weight: "<<Nmoca_reco_sel_no_weight<<endl;
   if(!isData) cout<<"selected events in mc at reco level - only vertex weight: "<<Nmoca_reco_sel_no_pu_weight<<endl;
   if(!isData) cout<<"selected events in mc at reco level - only pu weight: "<<Nmoca_reco_sel_no_vz_weight<<endl<<endl;
   
   long double ratio_moca_reco_sel_1 = 0;						       
   if(!isData) ratio_moca_reco_sel_1 = Nmoca_reco_sel_no_weight/Nmoca_reco_sel;		   	       
   if(!isData) cout<<"no weight / global weight: "<<ratio_moca_reco_sel_1<<endl;              
   
   long double ratio_moca_reco_sel_2 = 0;						       
   if(!isData) ratio_moca_reco_sel_2 = Nmoca_reco_sel_no_weight/Nmoca_reco_sel_no_pu_weight;  	       
   if(!isData) cout<<"no weight / only vertex weight: "<<ratio_moca_reco_sel_2<<endl;         
   
   long double ratio_moca_reco_sel_3 = 0;						       
   if(!isData) ratio_moca_reco_sel_3 = Nmoca_reco_sel_no_weight/Nmoca_reco_sel_no_vz_weight;  	       
   if(!isData) cout<<"no weight / only pu weight: "<<ratio_moca_reco_sel_3<<endl<<endl;       
   
   TString temp = outputfile_name;
   TString logfile_name(temp.Remove(temp.Length()-4));
   logfile_name+="txt";
   
   ofstream logfile;
   logfile.open(logfile_name);
   
   if(isData) {
     logfile<<endl<<"data selection: "<<endl;
     logfile<<"number of events before selection: "<<hselection_data->GetBinContent(1)<<endl;
     logfile<<"number of events after run and lumi section: "<<hselection_data->GetBinContent(2)<<endl;
     logfile<<"number of events after L1TT selection: "<<hselection_data->GetBinContent(3)<<endl;
     logfile<<"number of events after HLT: "<<hselection_data->GetBinContent(4)<<endl;
     logfile<<"number of events after at least one good vertex: "<<hselection_data->GetBinContent(5)<<endl;
     logfile<<"number of events after exactly one good vertex: "<<hselection_data->GetBinContent(6)<<endl;
   }
   
   if(!isData) {
     logfile<<endl<<"moca at reco level: "<<endl;
     logfile<<"number of events before selection: "<<hselection_moca_reco->GetBinContent(1)<<endl;
     logfile<<"number of events after at least one good vertex: "<<hselection_moca_reco->GetBinContent(2)<<endl;
     logfile<<"number of events after exactly one good vertex: "<<hselection_moca_reco->GetBinContent(3)<<endl;
   }
   
   logfile<<endl<<endl;
   logfile<<"run on "<<file_nb<<" file(s)"<<endl<<endl;
   
   logfile<<"total number of events - global weight: "<<Nevt_tot<<endl;
   logfile<<"total number of events - no weight: "<<Nevt_tot_no_weight<<endl;
   logfile<<"total number of events - only vertex weight: "<<Nevt_tot_no_pu_weight<<endl;
   logfile<<"total number of events - only pu weight: "<<Nevt_tot_no_vz_weight<<endl<<endl;
   
   logfile<<"no weight / global weight: "<<ratio_evt_tot_1<<endl;
   logfile<<"no weight / only vertex weight: "<<ratio_evt_tot_2<<endl;
   logfile<<"no weight / only pu weight: "<<ratio_evt_tot_3<<endl<<endl;
   
   if(isData) logfile<<"selected events in data - global weight: "<<Ndata_sel<<endl;
   if(isData) logfile<<"selected events in data - no weight: "<<Ndata_sel_no_weight<<endl;
   if(isData) logfile<<"selected events in data - only vertex weight: "<<Ndata_sel_no_pu_weight<<endl;
   if(isData) logfile<<"selected events in data - only pu weight: "<<Ndata_sel_no_vz_weight<<endl<<endl;
   
   if(isData) logfile<<"no weight / global weight: "<<ratio_data_sel_1<<endl;
   if(isData) logfile<<"no weight / only vertex weight: "<<ratio_data_sel_2<<endl;
   if(isData) logfile<<"no weight / only pu weight: "<<ratio_data_sel_3<<endl<<endl;
   
   if(!isData) logfile<<"selected events in mc at reco level - global weight: "<<Nmoca_reco_sel<<endl;
   if(!isData) logfile<<"selected events in mc at reco level - no weight: "<<Nmoca_reco_sel_no_weight<<endl;
   if(!isData) logfile<<"selected events in mc at reco level - only vertex weight: "<<Nmoca_reco_sel_no_pu_weight<<endl;
   if(!isData) logfile<<"selected events in mc at reco level - only pu weight: "<<Nmoca_reco_sel_no_vz_weight<<endl<<endl;
   
   if(!isData) logfile<<"no weight / global weight: "<<ratio_moca_reco_sel_1<<endl;
   if(!isData) logfile<<"no weight / only vertex weight: "<<ratio_moca_reco_sel_2<<endl;
   if(!isData) logfile<<"no weight / only pu weight: "<<ratio_moca_reco_sel_3<<endl<<endl;
   
   logfile.close();
   cout<<endl<<"log file: "<<logfile_name<<" created"<<endl;
   
   
   //-- if(!isData) cout<<Nmoca_reco_sel_HFplus10<<" selected events in mc at reco level with activity in HF plus 10 degrees"<<endl;
   //-- if(!isData) cout<<Nmoca_reco_sel_HFplus20<<" selected events in mc at reco level with activity in HF plus 20 degrees"<<endl<<endl;
   
   //-- if(!isData) cout<<Nmoca_reco_sel_HFminus10<<" selected events in mc at reco level with activity in HF minus 10 degrees"<<endl;
   //-- if(!isData) cout<<Nmoca_reco_sel_HFminus20<<" selected events in mc at reco level with activity in HF minus 20 degrees"<<endl<<endl;
   
   //-- cout<<endl<<endl;
   //-- cout<<"max layer information"<<endl<<endl;
   //-- for (int idet = 0; idet < detsize; idet++) 
   //--   cout<<DetName[idet]<<" has a maximum number of "<<MaxLayer[idet]<<" "<<SubdetName[idet]<<endl;
}



