#include "TrackCommon.h"

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

TrackCommon::TrackCommon() { }

TrackCommon::~TrackCommon() { }

void TrackCommon::GetBranches(TTree* tree,bool isData) {

  //-- Event id
  b_Run = tree->GetBranch("Run");
  b_LumiSection = tree->GetBranch("LumiSection");

  //-- L1TT
  if(isData) {
    b_L1TTname = tree->GetBranch("L1TTname");    
    b_L1TTbit = tree->GetBranch("L1TTbit");	    
    b_L1TTdecision = tree->GetBranch("L1TTdecision");
  }

  //-- HLT 
  if(isData) {
    b_HLTname = tree->GetBranch("HLTname");     	
    b_HLTindex = tree->GetBranch("HLTindex");    
    b_HLTprescale = tree->GetBranch("HLTprescale"); 	
    b_L1prescale = tree->GetBranch("L1prescale");  	
    b_HLTdecision = tree->GetBranch("HLTdecision");      
  }
  
  //-- beam spot
  b_BSx = tree->GetBranch("BSx");
  b_BSy = tree->GetBranch("BSy");
  b_BSz = tree->GetBranch("BSz");

  //-- vertex
  b_nVertex = tree->GetBranch("nVertex");						      
  
  b_Vx = tree->GetBranch("Vx");
  b_Vy = tree->GetBranch("Vy");
  b_Vz = tree->GetBranch("Vz");				     
  b_Vrho = tree->GetBranch("Vrho");
  
  b_VxError = tree->GetBranch("VxError");			     
  b_VyError = tree->GetBranch("VyError");			     
  b_VzError = tree->GetBranch("VzError");			     
  
  b_Vchi2 = tree->GetBranch("Vchi2");				     
  b_Vndof = tree->GetBranch("Vndof");				     
  b_Vchi2ndof = tree->GetBranch("Vchi2ndof");			     
  
  b_VisFake = tree->GetBranch("VisFake");				     
  b_VisValid = tree->GetBranch("VisValid");				     
  b_VisBest = tree->GetBranch("VisBest");				     
  
  //-- tracks associated to the vertex
  b_VnTrack = tree->GetBranch("VnTrack");				     
  b_VsumPt = tree->GetBranch("VsumPt");				     
  
  b_VTrackPt = tree->GetBranch("VTrackPt");	     
  b_VTrackEta = tree->GetBranch("VTrackEta");	     
  b_VTrackPhi = tree->GetBranch("VTrackPhi");	     
  
  b_VTrackPtError = tree->GetBranch("VTrackPtError"); 	     
  b_VTrackEtaError = tree->GetBranch("VTrackEtaError");       
  b_VTrackPhiError = tree->GetBranch("VTrackPhiError");       
  
  b_VTrackdxy = tree->GetBranch("VTrackdxy");    	     
  b_VTrackdxyError = tree->GetBranch("VTrackdxyError");	     
  							     
  b_VTrackdz = tree->GetBranch("VTrackdz");     	     
  b_VTrackdzError = tree->GetBranch("VTrackdzError");	     
  							     
  b_VTrackCharge = tree->GetBranch("VTrackCharge");     	     
  							     
  b_VTrackchi2 = tree->GetBranch("VTrackchi2");    	     
  b_VTrackndof = tree->GetBranch("VTrackndof");    	     
  b_VTrackchi2ndof = tree->GetBranch("VTrackchi2ndof");	     
  							     
  b_VTrackValidHit = tree->GetBranch("VTrackValidHit"); 	     
  b_VTrackLostHit = tree->GetBranch("VTrackLostHit");  	     
  	
  b_VTrackValidPixelHit = tree->GetBranch("VTrackValidPixelHit");
  b_VTrackValidStripHit = tree->GetBranch("VTrackValidStripHit"); 
                                                     
  b_VTrackWeight = tree->GetBranch("VTrackWeight");                     
  	     
  b_VTrackQuality = tree->GetBranch("VTrackQuality");     

  //-- all tracks
  b_nTrack = tree->GetBranch("nTrack");
  					                                               
  b_TrackPt = tree->GetBranch("TrackPt");
  b_TrackEta = tree->GetBranch("TrackEta");
  b_TrackPhi = tree->GetBranch("TrackPhi");
  					  						 
  b_TrackPtError = tree->GetBranch("TrackPtError");
  b_TrackEtaError = tree->GetBranch("TrackEtaError");
  b_TrackPhiError = tree->GetBranch("TrackPhiError");
  					  						 
  b_Trackdxy = tree->GetBranch("Trackdxy");
  b_TrackdxyError = tree->GetBranch("TrackdxyError");
  					  						 
  b_Trackdz = tree->GetBranch("Trackdz");
  b_TrackdzError = tree->GetBranch("TrackdzError");
  					  						 
  b_TrackCharge = tree->GetBranch("TrackCharge");
  					  						 
  b_Trackchi2 = tree->GetBranch("Trackchi2");
  b_Trackndof = tree->GetBranch("Trackndof");
  b_Trackchi2ndof = tree->GetBranch("Trackchi2ndof");
  					  						 
  b_TrackValidHit = tree->GetBranch("TrackValidHit");
  b_TrackLostHit = tree->GetBranch("TrackLostHit");
  					  						 
  b_TrackQuality = tree->GetBranch("TrackQuality");
  
  //-- all track rechits
  b_TrackHitR = tree->GetBranch("TrackHitR");      	
  b_TrackHitTheta = tree->GetBranch("TrackHitTheta");    
  b_TrackHitPhi = tree->GetBranch("TrackHitPhi");           
  							                                                             
  b_TrackHitEta = tree->GetBranch("TrackHitEta");           
  							                                                             
  b_TrackHitX = tree->GetBranch("TrackHitX");      	
  b_TrackHitY = tree->GetBranch("TrackHitY");      	
  b_TrackHitZ = tree->GetBranch("TrackHitZ");        
  							                                                             
  b_TrackHitDet = tree->GetBranch("TrackHitDet");         
  b_TrackHitSubdet = tree->GetBranch("TrackHitSubdet");   
  
  //-- pileup
  if(!isData) {
    b_nbxint = tree->GetBranch("nbxint");
    b_nbxearly = tree->GetBranch("nbxearly");
    b_nbxlate = tree->GetBranch("nbxlate");
    b_nbxtot = tree->GetBranch("nbxtot");	 
   
    b_PUint = tree->GetBranch("PUint");   
    b_PUearly = tree->GetBranch("PUearly"); 
    b_PUlate = tree->GetBranch("PUlate");  
    b_PUtot = tree->GetBranch("PUtot");   
			      
    b_NumInteraction = tree->GetBranch("NumInteraction");
  }             
    
  //-- generated particles
  if(!isData) {
    b_ParticleEnergy = tree->GetBranch("ParticleEnergy");	
    b_ParticlePt = tree->GetBranch("ParticlePt");		
    b_ParticleEta = tree->GetBranch("ParticleEta");		
    b_ParticlePhi = tree->GetBranch("ParticlePhi");		
      
    b_ParticleCharge = tree->GetBranch("ParticleCharge");	
    b_ParticleStatus = tree->GetBranch("ParticleStatus");	
    b_ParticleId = tree->GetBranch("ParticleId");              
  }                                                          

  //-- HF tower                                                                                                                   
  b_nHFTower = tree->GetBranch("nHFTower");
  b_HFTowerEtot = tree->GetBranch("HFTowerEtot");

  b_HFTowerEnergy = tree->GetBranch("HFTowerEnergy");
  b_HFTowerEta = tree->GetBranch("HFTowerEta");
  b_HFTowerPhi = tree->GetBranch("HFTowerPhi");

  b_HFplusTowerLeadingEnergy = tree->GetBranch("HFplusTowerLeadingEnergy");   	 
  b_HFplusTowerLeadingEta = tree->GetBranch("HFplusTowerLeadingEta");      	 
  b_HFplusTowerLeadingPhi = tree->GetBranch("HFplusTowerLeadingPhi");        

  b_HFminusTowerLeadingEnergy = tree->GetBranch("HFminusTowerLeadingEnergy");  
  b_HFminusTowerLeadingEta = tree->GetBranch("HFminusTowerLeadingEta");      	
  b_HFminusTowerLeadingPhi = tree->GetBranch("HFminusTowerLeadingPhi");        

  //-- HF rechit
  b_nHFRechit = tree->GetBranch("nHFRechit");	     
  	                        
  b_HFRechitEnergy = tree->GetBranch("HFRechitEnergy");  
  b_HFRechitEt = tree->GetBranch("HFRechitEt");	     
  b_HFRechitTime = tree->GetBranch("HFRechitTime");    
  	                        
  b_HFRechitiEta = tree->GetBranch("HFRechitiEta");    
  b_HFRechitiPhi = tree->GetBranch("HFRechitiPhi");    
  b_HFRechitDepth = tree->GetBranch("HFRechitDepth");   
  
  b_HFRechitEta = tree->GetBranch("HFRechitEta");     
  b_HFRechitPhi = tree->GetBranch("HFRechitPhi");  
}

void TrackCommon::SetBranchAddresses(TTree* tree,bool isData) {

  //-- Event id
  b_Run->SetAddress(&Run);
  b_LumiSection->SetAddress(&LumiSection);

  //-- L1TT
  if(isData) {
    b_L1TTname->SetAddress(&L1TTname);    
    b_L1TTbit->SetAddress(&L1TTbit);	    
    b_L1TTdecision->SetAddress(&L1TTdecision);
  }
  
  //-- HLT
  if(isData) {
    b_HLTname->SetAddress(&HLTname);     	  
    b_HLTindex->SetAddress(&HLTindex);    	  
    b_HLTprescale->SetAddress(&HLTprescale); 	  
    b_L1prescale->SetAddress(&L1prescale);  	  
    b_HLTdecision->SetAddress(&HLTdecision);      
  }
  
  //-- beam spot
  b_BSx->SetAddress(&BSx);    
  b_BSy->SetAddress(&BSy);
  b_BSz->SetAddress(&BSz);       

  //-- vertex
  b_nVertex->SetAddress(&nVertex);		     
    
  b_Vx->SetAddress(&Vx);			     
  b_Vy->SetAddress(&Vy);			     
  b_Vz->SetAddress(&Vz);			     
  b_Vrho->SetAddress(&Vrho);
  
  b_VxError->SetAddress(&VxError);		     
  b_VyError->SetAddress(&VyError);		     
  b_VzError->SetAddress(&VzError);		     
    
  b_Vchi2->SetAddress(&Vchi2);		     
  b_Vndof->SetAddress(&Vndof);		     
  b_Vchi2ndof->SetAddress(&Vchi2ndof);	     
  
  b_VisFake->SetAddress(&VisFake);		     
  b_VisValid->SetAddress(&VisValid);		     
  b_VisBest->SetAddress(&VisBest);	

  //-- tracks associated to the vertex
  b_VnTrack->SetAddress(&VnTrack);		     
  b_VsumPt->SetAddress(&VsumPt);		     
    
  b_VTrackPt->SetAddress(&VTrackPt);	     	     
  b_VTrackEta->SetAddress(&VTrackEta);	     
  b_VTrackPhi->SetAddress(&VTrackPhi);	     
    
  b_VTrackPtError->SetAddress(&VTrackPtError);     
  b_VTrackEtaError->SetAddress(&VTrackEtaError);   
  b_VTrackPhiError->SetAddress(&VTrackPhiError);   
  
  b_VTrackdxy->SetAddress(&VTrackdxy);    	     
  b_VTrackdxyError->SetAddress(&VTrackdxyError);   
    
  b_VTrackdz->SetAddress(&VTrackdz);     	     
  b_VTrackdzError->SetAddress(&VTrackdzError);     
    
  b_VTrackCharge->SetAddress(&VTrackCharge);       
    
  b_VTrackchi2->SetAddress(&VTrackchi2);    	     
  b_VTrackndof->SetAddress(&VTrackndof);    	     
  b_VTrackchi2ndof->SetAddress(&VTrackchi2ndof);   
    
  b_VTrackValidHit->SetAddress(&VTrackValidHit);   
  b_VTrackLostHit->SetAddress(&VTrackLostHit);     
    
  b_VTrackValidPixelHit->SetAddress(&VTrackValidPixelHit);	   
  b_VTrackValidStripHit->SetAddress(&VTrackValidStripHit); 	       
    
  b_VTrackWeight->SetAddress(&VTrackWeight);                         
 
  b_VTrackQuality->SetAddress(&VTrackQuality);    

  //-- all tracks
  b_nTrack->SetAddress(&nTrack);
    					                                               
  b_TrackPt->SetAddress(&TrackPt);
  b_TrackEta->SetAddress(&TrackEta);
  b_TrackPhi->SetAddress(&TrackPhi);
    					  						 
  b_TrackPtError->SetAddress(&TrackPtError);
  b_TrackEtaError->SetAddress(&TrackEtaError);
  b_TrackPhiError->SetAddress(&TrackPhiError);
  
  b_Trackdxy->SetAddress(&Trackdxy);
  b_TrackdxyError->SetAddress(&TrackdxyError);
  
  b_Trackdz->SetAddress(&Trackdz);
  b_TrackdzError->SetAddress(&TrackdzError);
    					  						 
  b_TrackCharge->SetAddress(&TrackCharge);
  
  b_Trackchi2->SetAddress(&Trackchi2);
  b_Trackndof->SetAddress(&Trackndof);
  b_Trackchi2ndof->SetAddress(&Trackchi2ndof);
  
  b_TrackValidHit->SetAddress(&TrackValidHit);
  b_TrackLostHit->SetAddress(&TrackLostHit);
    					  						 
  b_TrackQuality->SetAddress(&TrackQuality);

  //-- all track rechits
  b_TrackHitR->SetAddress(&TrackHitR);      	
  b_TrackHitTheta->SetAddress(&TrackHitTheta);    
  b_TrackHitPhi->SetAddress(&TrackHitPhi);           
    							                                                             
  b_TrackHitEta->SetAddress(&TrackHitEta);           
    							                                                             
  b_TrackHitX->SetAddress(&TrackHitX);      	
  b_TrackHitY->SetAddress(&TrackHitY);      	
  b_TrackHitZ->SetAddress(&TrackHitZ);        
  
  b_TrackHitDet->SetAddress(&TrackHitDet);         
  b_TrackHitSubdet->SetAddress(&TrackHitSubdet);      
  
  //-- pileup
  if(!isData) {
    b_nbxint->SetAddress(&nbxint);	 		    
    b_nbxearly->SetAddress(&nbxearly);	 		    
    b_nbxlate->SetAddress(&nbxlate);	 		    
    b_nbxtot->SetAddress(&nbxtot);	 		    
    
    b_PUint->SetAddress(&PUint);   		    
    b_PUearly->SetAddress(&PUearly); 		    
    b_PUlate->SetAddress(&PUlate);  		    
    b_PUtot->SetAddress(&PUtot);   		    
    
    b_NumInteraction->SetAddress(&NumInteraction); 
  }                                                       
  
  //-- generated particles
  if(!isData) {
    b_ParticleEnergy->SetAddress(&ParticleEnergy);	
    b_ParticlePt->SetAddress(&ParticlePt);		
    b_ParticleEta->SetAddress(&ParticleEta);		
    b_ParticlePhi->SetAddress(&ParticlePhi);		
      
    b_ParticleCharge->SetAddress(&ParticleCharge);	
    b_ParticleStatus->SetAddress(&ParticleStatus);	
    b_ParticleId->SetAddress(&ParticleId);              
  }
   
  //-- HF tower
  b_nHFTower->SetAddress(&nHFTower);
  b_HFTowerEtot->SetAddress(&HFTowerEtot);

  b_HFTowerEnergy->SetAddress(&HFTowerEnergy);
  b_HFTowerEta->SetAddress(&HFTowerEta);
  b_HFTowerPhi->SetAddress(&HFTowerPhi);       

  b_HFplusTowerLeadingEnergy->SetAddress(&HFplusTowerLeadingEnergy); 
  b_HFplusTowerLeadingEta->SetAddress(&HFplusTowerLeadingEta);       
  b_HFplusTowerLeadingPhi->SetAddress(&HFplusTowerLeadingPhi);       

  b_HFminusTowerLeadingEnergy->SetAddress(&HFminusTowerLeadingEnergy); 
  b_HFminusTowerLeadingEta->SetAddress(&HFminusTowerLeadingEta);       
  b_HFminusTowerLeadingPhi->SetAddress(&HFminusTowerLeadingPhi);       

  //-- HF rechit
  b_nHFRechit->SetAddress(&nHFRechit);	     
    	                        
  b_HFRechitEnergy->SetAddress(&HFRechitEnergy);  
  b_HFRechitEt->SetAddress(&HFRechitEt);	     
  b_HFRechitTime->SetAddress(&HFRechitTime);    
    	                        
  b_HFRechitiEta->SetAddress(&HFRechitiEta);    
  b_HFRechitiPhi->SetAddress(&HFRechitiPhi);    
  b_HFRechitDepth->SetAddress(&HFRechitDepth);   
    	                        
  b_HFRechitEta->SetAddress(&HFRechitEta);     
  b_HFRechitPhi->SetAddress(&HFRechitPhi);          
}

void TrackCommon::GetTreeEntries(int ievt, bool isData) {

  //-- Event id
  b_Run->GetEntry(ievt);    
  b_LumiSection->GetEntry(ievt);    

  //-- L1TT
  if(isData) {
    b_L1TTname->GetEntry(ievt);    
    b_L1TTbit->GetEntry(ievt);	    
    b_L1TTdecision->GetEntry(ievt);
  }

  //-- HLT
  if(isData) {
    b_HLTname->GetEntry(ievt);    
    b_HLTindex->GetEntry(ievt);   
    b_HLTprescale->GetEntry(ievt);
    b_L1prescale->GetEntry(ievt); 
    b_HLTdecision->GetEntry(ievt);
  }

  //-- beam spot
  b_BSx->GetEntry(ievt);       
  b_BSy->GetEntry(ievt);	   
  b_BSz->GetEntry(ievt);       

  //-- vertex
  b_nVertex->GetEntry(ievt);
    
  b_Vx->GetEntry(ievt);
  b_Vy->GetEntry(ievt);
  b_Vz->GetEntry(ievt);
  b_Vrho->GetEntry(ievt);
  
  b_VxError->GetEntry(ievt);
  b_VyError->GetEntry(ievt);
  b_VzError->GetEntry(ievt);
      
  b_Vchi2->GetEntry(ievt);
  b_Vndof->GetEntry(ievt);
  b_Vchi2ndof->GetEntry(ievt);
  
  b_VisFake->GetEntry(ievt);
  b_VisValid->GetEntry(ievt);
  b_VisBest->GetEntry(ievt);
  
  //-- tracks associated to the vertex
  b_VnTrack->GetEntry(ievt);
  b_VsumPt->GetEntry(ievt);
  
  b_VTrackPt->GetEntry(ievt);
  b_VTrackEta->GetEntry(ievt);
  b_VTrackPhi->GetEntry(ievt);
  
  b_VTrackPtError->GetEntry(ievt);
  b_VTrackEtaError->GetEntry(ievt);
  b_VTrackPhiError->GetEntry(ievt);
  
  b_VTrackdxy->GetEntry(ievt);
  b_VTrackdxyError->GetEntry(ievt);
      
  b_VTrackdz->GetEntry(ievt);
  b_VTrackdzError->GetEntry(ievt);
  
  b_VTrackCharge->GetEntry(ievt);
  
  b_VTrackchi2->GetEntry(ievt);
  b_VTrackndof->GetEntry(ievt);
  b_VTrackchi2ndof->GetEntry(ievt);
  
  b_VTrackValidHit->GetEntry(ievt);
  b_VTrackLostHit->GetEntry(ievt);
  
  b_VTrackValidPixelHit->GetEntry(ievt);
  b_VTrackValidStripHit->GetEntry(ievt);
  
  b_VTrackWeight->GetEntry(ievt);

  b_VTrackQuality->GetEntry(ievt);
      
  //-- all tracks
  b_nTrack->GetEntry(ievt);
      
  b_TrackPt->GetEntry(ievt);
  b_TrackEta->GetEntry(ievt);
  b_TrackPhi->GetEntry(ievt);
      
  b_TrackPtError->GetEntry(ievt);
  b_TrackEtaError->GetEntry(ievt);
  b_TrackPhiError->GetEntry(ievt);
  
  b_Trackdxy->GetEntry(ievt);
  b_TrackdxyError->GetEntry(ievt);
      
  b_Trackdz->GetEntry(ievt);
  b_TrackdzError->GetEntry(ievt);
      
  b_TrackCharge->GetEntry(ievt);
      
  b_Trackchi2->GetEntry(ievt);
  b_Trackndof->GetEntry(ievt);
  b_Trackchi2ndof->GetEntry(ievt);
      
  b_TrackValidHit->GetEntry(ievt);
  b_TrackLostHit->GetEntry(ievt);
      
  b_TrackQuality->GetEntry(ievt);
      
  //-- all track rechits 
  b_TrackHitR->GetEntry(ievt);
  b_TrackHitTheta->GetEntry(ievt);
  b_TrackHitPhi->GetEntry(ievt);
    							                                                             
  b_TrackHitEta->GetEntry(ievt);
      
  b_TrackHitX->GetEntry(ievt);
  b_TrackHitY->GetEntry(ievt);
  b_TrackHitZ->GetEntry(ievt);
  
  b_TrackHitDet->GetEntry(ievt);
  b_TrackHitSubdet->GetEntry(ievt);

  //-- pileup
  if(!isData) {
    b_nbxint->GetEntry(ievt);	 		     
    b_nbxearly->GetEntry(ievt);	 		     				     
    b_nbxlate->GetEntry(ievt);	 		     
    b_nbxtot->GetEntry(ievt);	 		     
    
    b_PUint->GetEntry(ievt);
    b_PUearly->GetEntry(ievt);
    b_PUlate->GetEntry(ievt);
    b_PUtot->GetEntry(ievt);
    
    b_NumInteraction->GetEntry(ievt);
  }                                                
  
  //-- generated particles
  if(!isData) {
    b_ParticleEnergy->GetEntry(ievt);
    b_ParticlePt->GetEntry(ievt);
    b_ParticleEta->GetEntry(ievt);
    b_ParticlePhi->GetEntry(ievt);
    
    b_ParticleCharge->GetEntry(ievt);
    b_ParticleStatus->GetEntry(ievt);
    b_ParticleId->GetEntry(ievt);
  }

  //-- HF tower                                
  b_nHFTower->GetEntry(ievt);
  b_HFTowerEtot->GetEntry(ievt);
  
  b_HFTowerEnergy->GetEntry(ievt);
  b_HFTowerEta->GetEntry(ievt);
  b_HFTowerPhi->GetEntry(ievt);
  
  b_HFplusTowerLeadingEnergy->GetEntry(ievt);
  b_HFplusTowerLeadingEta->GetEntry(ievt);
  b_HFplusTowerLeadingPhi->GetEntry(ievt);    
  
  b_HFminusTowerLeadingEnergy->GetEntry(ievt); 
  b_HFminusTowerLeadingEta->GetEntry(ievt);	  
  b_HFminusTowerLeadingPhi->GetEntry(ievt);    
      
  //-- HF rechit						      
  b_nHFRechit->GetEntry(ievt);
  
  b_HFRechitEnergy->GetEntry(ievt); 
  b_HFRechitEt->GetEntry(ievt);
  b_HFRechitTime->GetEntry(ievt);   
  
  b_HFRechitiEta->GetEntry(ievt);   
  b_HFRechitiPhi->GetEntry(ievt);   
  b_HFRechitDepth->GetEntry(ievt);  
  
  b_HFRechitEta->GetEntry(ievt);    
  b_HFRechitPhi->GetEntry(ievt);
}


   

    
