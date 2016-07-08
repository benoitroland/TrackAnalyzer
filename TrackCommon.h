#ifndef TrackCommon_h
#define TrackCommon_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TBranch.h>

#include "Common.h"

class TrackCommon: public Common {

 public:

  TrackCommon();
  virtual ~TrackCommon();

  void GetBranches(TTree* tree,bool isData);
  void SetBranchAddresses(TTree* tree,bool isData);
  void GetTreeEntries(int ievt, bool isData);

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

  //-- all tracks
  int nTrack = 0;                        
  
  std::vector<double> *TrackPt = NULL;                 
  std::vector<double> *TrackEta = NULL;                
  std::vector<double> *TrackPhi = NULL;                
  
  std::vector<double> *TrackPtError = NULL;            
  std::vector<double> *TrackEtaError = NULL;           
  std::vector<double> *TrackPhiError = NULL;           
  
  std::vector<double> *Trackdxy = NULL; 
  std::vector<double> *TrackdxyError = NULL;
  
  std::vector<double> *Trackdz = NULL;
  std::vector<double> *TrackdzError = NULL;
  
  std::vector<int> *TrackCharge = NULL;                
  
  std::vector<double> *Trackchi2 = NULL;
  std::vector<double> *Trackndof = NULL;
  std::vector<double> *Trackchi2ndof = NULL;
  
  std::vector<int> *TrackValidHit = NULL;
  std::vector<int> *TrackLostHit = NULL;
  
  std::vector<int> *TrackQuality = NULL;    

  //-- all track rechits
  std::vector< std::vector<double> > *TrackHitR = NULL;      
  std::vector< std::vector<double> > *TrackHitTheta = NULL;     
  std::vector< std::vector<double> > *TrackHitPhi = NULL;            
  
  std::vector< std::vector<double> > *TrackHitEta = NULL;            
  
  std::vector< std::vector<double> > *TrackHitX = NULL;      
  std::vector< std::vector<double> > *TrackHitY = NULL;      
  std::vector< std::vector<double> > *TrackHitZ = NULL;         
  
  std::vector< std::vector<int> > *TrackHitDet = NULL;          
  std::vector< std::vector<int> > *TrackHitSubdet = NULL;  

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
  
  double HFplusTowerLeadingEnergy = 0;	 
  double HFplusTowerLeadingEta = 0;	 
  double HFplusTowerLeadingPhi = 0;      
  
  double HFminusTowerLeadingEnergy = 0;   
  double HFminusTowerLeadingEta = 0;	   
  double HFminusTowerLeadingPhi = 0;      
  
  //-- HF rechit
  int nHFRechit = 0;
  
  std::vector<double> *HFRechitEnergy = NULL;
  std::vector<double> *HFRechitEt = NULL;
  std::vector<double> *HFRechitTime = NULL;
  
  std::vector<int> *HFRechitiEta = NULL;
  std::vector<int> *HFRechitiPhi = NULL;
  std::vector<int> *HFRechitDepth = NULL;
  
  std::vector<double> *HFRechitEta = NULL;
  std::vector<double> *HFRechitPhi = NULL;     

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

  //-- all tracks
  TBranch *b_nTrack = NULL; 
  						                                               
  TBranch *b_TrackPt = NULL;  
  TBranch *b_TrackEta = NULL; 
  TBranch *b_TrackPhi = NULL; 
  
  TBranch *b_TrackPtError = NULL;  
  TBranch *b_TrackEtaError = NULL; 
  TBranch *b_TrackPhiError = NULL; 
    						  						 
  TBranch *b_Trackdxy = NULL;      
  TBranch *b_TrackdxyError = NULL; 
  
  TBranch *b_Trackdz = NULL; 
  TBranch *b_TrackdzError = NULL; 
    						  						 
  TBranch *b_TrackCharge = NULL; 
    						  						 
  TBranch *b_Trackchi2 = NULL; 
  TBranch *b_Trackndof = NULL; 
  TBranch *b_Trackchi2ndof = NULL; 
    						  						 
  TBranch *b_TrackValidHit = NULL; 
  TBranch *b_TrackLostHit = NULL; 
    						  						 
  TBranch *b_TrackQuality = NULL; 

  //-- all track rechits
  TBranch *b_TrackHitR = NULL;     
  TBranch *b_TrackHitTheta = NULL; 
  TBranch *b_TrackHitPhi = NULL; 
  								                                                             
  TBranch *b_TrackHitEta = NULL; 
  								                                                             
  TBranch *b_TrackHitX = NULL; 	
  TBranch *b_TrackHitY = NULL; 	
  TBranch *b_TrackHitZ = NULL; 
  								                                                             
  TBranch *b_TrackHitDet = NULL; 
  TBranch *b_TrackHitSubdet = NULL; 

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

  //-- HF rechit
  TBranch *b_nHFRechit = NULL; 
  
  TBranch *b_HFRechitEnergy = NULL; 
  TBranch *b_HFRechitEt = NULL;        
  TBranch *b_HFRechitTime = NULL;   
  
  TBranch *b_HFRechitiEta = NULL;   
  TBranch *b_HFRechitiPhi = NULL;   
  TBranch *b_HFRechitDepth = NULL;  
  
  TBranch *b_HFRechitEta = NULL;   
  TBranch *b_HFRechitPhi = NULL; 

 private:
};

#endif

