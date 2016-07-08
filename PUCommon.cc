#include "PUCommon.h"

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

PUCommon::PUCommon() { }

PUCommon::~PUCommon() { }

void PUCommon::GetBranches(TTree* tree,bool isData) {

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
    
}

void PUCommon::SetBranchAddresses(TTree* tree,bool isData) {

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
  
}

void PUCommon::GetTreeEntries(int ievt, bool isData) {

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
  
}


  







  
