#ifndef RhoCommon_h
#define RhoCommon_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TTree.h>
#include <TBranch.h>

#include "Common.h"

class RhoCommon: public Common {

 public:

  RhoCommon();
  virtual ~RhoCommon();

  void GetBranches(TTree* tree,bool isData);
  void SetBranchAddresses(TTree* tree,bool isData);
  void GetTreeEntries(int ievt, bool isData);

  //-- constant

  
  static const int nHFcut = 3;
  static const int nHFsel = 4;
    
  static const int neta_edge = 5;
  const double delta_eta_0 = 0.4;   
  static const int n_delta_eta = 11;

  TString HFcut[3] = {"cut down","cut nominal","cut up"};
  TString HFsel[4] = {"inclusive selection","inelastic selection","NSD selection","SD selection"};

  //-- HF eta range                                                                                                                         
  bool GetHFrange(double eta);

  //-- HF activity                                                                                                                                      
  bool GetHFActivity(int status, double energy, double eta);

  //-- good track decision                                                                                                                       
  bool GetGoodTrack(int quality,double pTerror,double Sigxy,double Sigz,int npixelhit,double eta, double phi);

  //-- good vertex decision                                                                                                                       
  bool GetGoodVertex15(bool isfake,bool isvalid,double vz,double rho,int nGoodTrack);

  //-- good vertex decision                                                                                                            
  bool GetGoodVertex20(bool isfake,bool isvalid,double vz,double rho,int nGoodTrack);

  //-- good particle decision
  bool GetGoodParticle(int status, int charge, double eta, double pT);

  //-- compute correlation coefficient                                                                                      
  void GetCorrelation(TH1D* hnFwd[n_delta_eta], TH1D* hnBwd[n_delta_eta], TH1D* hnFwdBwd[n_delta_eta], TH1D* hrho);
  void GetCorrelation(TH2D* hcorrelation[n_delta_eta], TH1D* hrho);
  void GetCorrelation(TH2D* hcorrelation[nHFsel][n_delta_eta], TH1D* hrho[nHFsel]);

  //-- compute histo ratio
  void ComputeHistoRatio(TH1D* hup[nHFsel], TH1D* hdown[nHFsel], TH1D* hratio[nHFsel]);

  //-- create histograms
  void CreateHistoSelection();
  
  void CreateHistoHFTower();
  void CreateHistoHFplusTower();
  void CreateHistoHFminusTower();

  void CreateHistoTrack();
  void CreateHistoVertex();
  void CreateHistoBS();
  
  void CreateHistoPU();
  void CreateHistoRho();

  //-- check histograms
  void CheckHistoSelection();

  void CheckHistoHFTower();
  void CheckHistoHFplusTower();
  void CheckHistoHFminusTower();

  void CheckHistoTrack();
  void CheckHistoVertex();
  void CheckHistoBS();

  void CheckHistoPU(bool isData);
  void CheckHistoRho(bool isData);
  
  //-- write histograms
  void WriteHistoSelection();

  void WriteHistoHFTower();
  void WriteHistoHFplusTower();
  void WriteHistoHFminusTower();

  void WriteHistoTrack();
  void WriteHistoVertex();
  void WriteHistoBS();

  void WriteHistoPU(bool isData);
  void WriteHistoRho(bool isData);

  //---------------//
  //-- variables --//
  //---------------//

  //-- Event id
  int Run = 0;
  int LumiSection = 0;

  //-- L1TT
  std::vector<std::string> *L1TTname = NULL;
  std::vector<int> *L1TTbit = NULL;
  std::vector<bool> *L1TTdecision = NULL;

  //-- HLT
  std::vector<std::string> *HLTname = NULL;
  std::vector<unsigned int> *HLTindex = NULL;
  std::vector<double> *HLTprescale = NULL;
  std::vector<double> *L1prescale = NULL;
  std::vector<bool> *HLTdecision = NULL;      

  //-- beam spot
  double BSx = 0;
  double BSy = 0;
  double BSz = 0;

  //-- vertex 
  int nVertex = 0;
  
  std::vector<double> *Vx = NULL;
  std::vector<double> *Vy = NULL;
  std::vector<double> *Vz = NULL;
  std::vector<double> *Vrho = NULL;
  
  std::vector<double> *VxError = NULL;
  std::vector<double> *VyError = NULL;
  std::vector<double> *VzError = NULL;
  
  std::vector<double> *Vchi2 = NULL;
  std::vector<double> *Vndof = NULL;
  std::vector<double> *Vchi2ndof = NULL;
  
  std::vector<bool> *VisFake = NULL;
  std::vector<bool> *VisValid = NULL;
  std::vector<bool> *VisBest = NULL;
  
  //-- tracks associated to the vertex
  std::vector<int> *VnTrack = NULL;
  std::vector<double> *VsumPt = NULL;
  
  std::vector< std::vector<double> > *VTrackPt = NULL;
  std::vector< std::vector<double> > *VTrackEta = NULL;
  std::vector< std::vector<double> > *VTrackPhi = NULL;
  
  std::vector< std::vector<double> > *VTrackPtError = NULL; 
  std::vector< std::vector<double> > *VTrackEtaError = NULL;  
  std::vector< std::vector<double> > *VTrackPhiError = NULL;  
  
  std::vector< std::vector<double> > *VTrackdxy = NULL;    
  std::vector< std::vector<double> > *VTrackdxyError = NULL;
  
  std::vector< std::vector<double> > *VTrackdz = NULL;     
  std::vector< std::vector<double> > *VTrackdzError = NULL;
  
  std::vector< std::vector<int> > *VTrackCharge = NULL;     
  
  std::vector< std::vector<double> > *VTrackchi2 = NULL;    
  std::vector< std::vector<double> > *VTrackndof = NULL;    
  std::vector< std::vector<double> > *VTrackchi2ndof = NULL;
  
  std::vector< std::vector<int> > *VTrackValidHit = NULL; 
  std::vector< std::vector<int> > *VTrackLostHit = NULL;  
  
  std::vector< std::vector<int> > *VTrackValidPixelHit = NULL;
  std::vector< std::vector<int> > *VTrackValidStripHit = NULL;
  
  std::vector< std::vector<double> > *VTrackWeight = NULL;                     
  
  std::vector< std::vector<int> > *VTrackQuality = NULL;     

  //-- pileup                                                                                                      
  int nbxint = 0;
  int nbxearly = 0;
  int nbxlate = 0;
  int nbxtot = 0;
    
  double PUint = 0;
  double PUearly = 0;
  double PUlate = 0;
  double PUtot = 0;         
  
  double NumInteraction = 0;
 
  //-- generated particles
  std::vector<double> *ParticleEnergy = NULL;
  std::vector<double> *ParticlePt = NULL;
  std::vector<double> *ParticleEta = NULL;
  std::vector<double> *ParticlePhi = NULL;

  std::vector<int> *ParticleCharge = NULL;
  std::vector<int> *ParticleStatus = NULL;
  std::vector<int> *ParticleId = NULL;

  //-- HF tower                                                                                                                                
  int nHFTower = 0;
  double HFTowerEtot = 0;
  
  std::vector<double> *HFTowerEnergy = NULL;
  std::vector<double> *HFTowerEta = NULL;
  std::vector<double> *HFTowerPhi = NULL;      
  
  //--------------//
  //-- branches --//
  //--------------//

  //-- Event id
  TBranch *b_Run = NULL;
  TBranch *b_LumiSection = NULL;

  //-- L1TT
  TBranch *b_L1TTname = NULL;
  TBranch *b_L1TTbit = NULL;
  TBranch *b_L1TTdecision = NULL;

  //-- HLT
  TBranch *b_HLTname = NULL;
  TBranch *b_HLTindex = NULL;
  TBranch *b_HLTprescale = NULL;
  TBranch *b_L1prescale = NULL;
  TBranch *b_HLTdecision = NULL;

  //-- beam spot
  TBranch *b_BSx = NULL;
  TBranch *b_BSy = NULL;
  TBranch *b_BSz = NULL;

  //-- vertex
  TBranch *b_nVertex = NULL;
      								     
  TBranch *b_Vx = NULL;
  TBranch *b_Vy = NULL;
  TBranch *b_Vz = NULL;
  TBranch *b_Vrho = NULL;

  TBranch *b_VxError = NULL; 
  TBranch *b_VyError = NULL; 
  TBranch *b_VzError = NULL; 
    								     
  TBranch *b_Vchi2 = NULL; 
  TBranch *b_Vndof = NULL; 
  TBranch *b_Vchi2ndof = NULL; 
    								     
  TBranch *b_VisFake = NULL; 
  TBranch *b_VisValid = NULL; 
  TBranch *b_VisBest = NULL;
    	
  //-- tracks associated to the vertex							     
  TBranch *b_VnTrack = NULL;
  TBranch *b_VsumPt = NULL; 
    								     
  TBranch *b_VTrackPt = NULL; 
  TBranch *b_VTrackEta = NULL;
  TBranch *b_VTrackPhi = NULL;
    								     
  TBranch *b_VTrackPtError = NULL; 
  TBranch *b_VTrackEtaError = NULL;
  TBranch *b_VTrackPhiError = NULL;
    								     
  TBranch *b_VTrackdxy = NULL; 
  TBranch *b_VTrackdxyError = NULL; 
    								     
  TBranch *b_VTrackdz = NULL; 
  TBranch *b_VTrackdzError = NULL; 
    								     
  TBranch *b_VTrackCharge = NULL;
    								     
  TBranch *b_VTrackchi2 = NULL; 
  TBranch *b_VTrackndof = NULL; 
  TBranch *b_VTrackchi2ndof = NULL; 
    								     
  TBranch *b_VTrackValidHit = NULL; 
  TBranch *b_VTrackLostHit = NULL; 
    		
  TBranch *b_VTrackValidPixelHit = NULL; 
  TBranch *b_VTrackValidStripHit = NULL; 
                                                              
  TBranch *b_VTrackWeight = NULL; 
  
  TBranch *b_VTrackQuality = NULL; 

  //-- pileup
  TBranch *b_nbxint = NULL;
  TBranch *b_nbxearly = NULL;
  TBranch *b_nbxlate = NULL;
  TBranch *b_nbxtot = NULL;	      
			      
  TBranch *b_PUint = NULL;	      
  TBranch *b_PUearly = NULL;	      
  TBranch *b_PUlate = NULL;	      
  TBranch *b_PUtot = NULL;	      
    		      
  TBranch *b_NumInteraction = NULL;

  //-- generated particles
  TBranch *b_ParticleEnergy = NULL;
  TBranch *b_ParticlePt = NULL;
  TBranch *b_ParticleEta = NULL;
  TBranch *b_ParticlePhi = NULL;

  TBranch *b_ParticleCharge = NULL;
  TBranch *b_ParticleStatus = NULL;
  TBranch *b_ParticleId = NULL;

  //-- HF tower                                                                                                                   
  TBranch *b_nHFTower = NULL; 
  TBranch *b_HFTowerEtot = NULL; 

  TBranch *b_HFTowerEnergy = NULL; 
  TBranch *b_HFTowerEta = NULL; 
  TBranch *b_HFTowerPhi = NULL; 

  TBranch *b_HFplusTowerLeadingEnergy = NULL; 
  TBranch *b_HFplusTowerLeadingEta = NULL;    
  TBranch *b_HFplusTowerLeadingPhi = NULL;    

  TBranch *b_HFminusTowerLeadingEnergy = NULL; 
  TBranch *b_HFminusTowerLeadingEta = NULL;    	
  TBranch *b_HFminusTowerLeadingPhi = NULL;    

  //----------------//
  //-- histograms --//
  //----------------//

  //-- histo label and title
  TString label,title,title_fwd,title_bwd;
  TString xleg, yleg;

  //-- selection
  TH1D *hselection_data = NULL;
  TH1D *hselection_moca_reco = NULL;
  TH1D *hselection_moca_gen = NULL;

  //-- HF energy cut
  double EHFcut[nHFcut] = {4,5,6};
  double EGencut = 5;

  //-- HF tower
  TH1D *hnHFTower[nHFcut];
  TH1D *hHFTowerEtot[nHFcut];
  TH1D *hHFTowerLeadingEnergy[nHFcut];

  TH1D *hHFTowerEnergy[nHFcut];
  TH1D *hHFTowerEta[nHFcut];
  TH1D *hHFTowerPhi[nHFcut];            

  TH2D *hHFTowerEtaPhi[nHFcut];
  
  //-- HF plus tower   
  TH1D *hnHFplusTower[nHFcut];   
  TH1D *hHFplusTowerEtot[nHFcut];   
  TH1D *hHFplusTowerLeadingEnergy[nHFcut];   
  
  TH1D *hHFplusTowerEnergy[nHFcut];   
  TH1D *hHFplusTowerEta[nHFcut];   
  TH1D *hHFplusTowerPhi[nHFcut];

  TH2D *hHFplusTowerEtaPhi[nHFcut];

  //-- HF minus tower
  TH1D *hnHFminusTower[nHFcut];   
  TH1D *hHFminusTowerEtot[nHFcut];   
  TH1D *hHFminusTowerLeadingEnergy[nHFcut];   
  
  TH1D *hHFminusTowerEnergy[nHFcut];   
  TH1D *hHFminusTowerEta[nHFcut];   
  TH1D *hHFminusTowerPhi[nHFcut];

  TH2D *hHFminusTowerEtaPhi[nHFcut];

  //-- track distributions exactly one good vertex
  TH1D *hTrackNum = NULL;
  TH1D *hTrackPt = NULL; 
  TH1D *hTrackEta = NULL;

  TH1D *hTrackPhi = NULL;
  TH1D *hTrackPhiCentral = NULL;
  TH1D *hTrackPhiFwd = NULL;

  TH2D* hTrackEtaPhi = NULL;

  TH1D *hTrackQuality = NULL;

  double eta_edge[neta_edge] = {0,0.6,1.2,1.8,2.4};
  
  TH1D *hTrackPtEta[neta_edge-1];
  
  TH1D *hTrackPtError = NULL;
  TH1D *hTrackEtaError = NULL;
  TH1D *hTrackPhiError = NULL;
  
  TProfile *hTrackPtErrorPt = NULL;
  TProfile *hTrackPtErrorEta = NULL;
  TProfile *hTrackPtErrorNhit = NULL;

  TH1D *hTrackdxy = NULL;
  TH1D *hTrackdz = NULL;

  TProfile *hTrackdxyEta = NULL;
  TProfile *hTrackdzEta = NULL;

  TProfile *hTrackdxyNhit = NULL;
  TProfile *hTrackdzNhit = NULL;
  
  TH2D *hTrackdxyEta2D = NULL;
  TH2D *hTrackdzEta2D = NULL;

  TH2D *hTrackdxyNum2D = NULL;
  TH2D *hTrackdzNum2D = NULL;

  TH1D *hTrackchi2 = NULL;
  TH1D *hTrackndof = NULL;
  TH1D *hTrackchi2ndof = NULL;

  TH1D *hTrackValidHit = NULL;
  TH1D *hTrackValidPixelHit = NULL;
  TH1D *hTrackValidStripHit = NULL;
 
  TH1D *hTrackWeight = NULL;

  TProfile *hTrackValidHitEta = NULL;
  TProfile *hTrackValidPixelHitEta = NULL;
  TProfile *hTrackValidStripHitEta = NULL;

  TH1D *hLeadingTrackPt = NULL;
  TH1D *hLeadingTrackEta = NULL;
  TH1D *hLeadingTrackPhi = NULL;

  //-- vertex distributions exactly one good vertex
  TH1D *hVx = NULL;
  TH1D *hVy = NULL;
  TH1D *hVz = NULL;
  TH2D *hVxy = NULL;
  TH1D *hVrho = NULL;

  TH1D *hVchi2 = NULL;
  TH1D *hVndof = NULL;
  TH1D *hVchi2ndof = NULL;

  TH1D *hVzError = NULL;
  TH2D *hVzErrorVz = NULL;
  TH2D *hVzErrorTrackNum = NULL;
  
  TH1D *hVNum = NULL;
  TH1D *hVNumGood = NULL;

  //-- beam spot distributions
  TH1D *hBSx = NULL;
  TH1D *hBSy = NULL;
  TH1D *hBSz = NULL;

  //-- pileup information
  TH1D *hnbxint = NULL;
  TH1D *hnbxearly = NULL;
  TH1D *hnbxlate = NULL;
  TH1D *hnbxtot = NULL;

  TH1D *hPUint = NULL;
  TH1D *hPUearly = NULL;
  TH1D *hPUlate = NULL;
  TH1D *hPUtot = NULL;

  TH1D* hNumInt = NULL;

  TH2D* hPUintVNum = NULL;
  TH2D* hPUearlyVNum = NULL;
  TH2D* hPUlateVNum = NULL;
  TH2D* hPUtotVNum = NULL;
  TH2D* hPUintVNum0 = NULL;

  //-- correlation coefficient  
  double eta_plus_min[n_delta_eta];
  double eta_plus_max[n_delta_eta];
  
  double eta_minus_min[n_delta_eta];
  double eta_minus_max[n_delta_eta];

  TH2D* hcorrelation[nHFsel][n_delta_eta];

  TH2D* hcorrelationDown[nHFsel][n_delta_eta];
  TH2D* hcorrelationUp[nHFsel][n_delta_eta];

  TH2D* hcorrelationPU[nHFsel][n_delta_eta];
  TH2D* hcorrelationEffi[nHFsel][n_delta_eta];

  TH2D* hcorrelationGen[nHFsel][n_delta_eta];
  TH2D* hcorrelationTheo[nHFsel][n_delta_eta];  

  TH1D* hrho[nHFsel];

  TH1D* hrhoDown[nHFsel];
  TH1D* hrhoUp[nHFsel];

  TH1D* hrhoPU[nHFsel];
  TH1D* hrhoEffi[nHFsel];

  TH1D* hrhoGen[nHFsel];
  TH1D* hrhoTheo[nHFsel];

  TH1D* hrhoCF[nHFsel];    //-- for nominal and Effi
  TH1D* hrhoCFPU[nHFsel];  //-- for PU

 private:
};

#endif

