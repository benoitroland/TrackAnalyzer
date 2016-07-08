#ifndef VertexCommon_h
#define VertexCommon_h

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

class VertexCommon: public Common {

 public:

  VertexCommon();
  virtual ~VertexCommon();

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

  //-------------//
  //-- branches -//
  //-------------//

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

 private:

};

#endif

