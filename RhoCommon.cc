#include "RhoCommon.h"

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

RhoCommon::RhoCommon() { }

RhoCommon::~RhoCommon() { }

  
//-- HF eta range
bool RhoCommon::GetHFrange(double eta) {

  bool ishf = true;

  if(TMath::Abs(eta) < 3) ishf = false;

  if(TMath::Abs(eta) > 5) ishf = false;

  return(ishf);
}

//-- HF activity
bool RhoCommon::GetHFActivity(int status, double energy, double eta) {

  bool activity = true;
			      
  if(status != 1) activity = false;
  if(energy < EGencut) activity = false;
  if(!GetHFrange(eta)) activity = false;

  return(activity);
}


//-- good track decision
bool RhoCommon::GetGoodTrack(int quality,double pTerror,double Sigxy,double Sigz,int npixelhit,double eta,double phi) {

  bool isgoodtrack = true;

  //-- quality = 2
  if(quality != 2) isgoodtrack = false;

  //-- pTerror <= 10
  if(pTerror > 10) isgoodtrack = false;
  
  //-- |Sigxy| <= 3 
  if(TMath::Abs(Sigxy) > 3) isgoodtrack = false;

  //-- |Sigz| <= 3
  if(TMath::Abs(Sigz) > 3) isgoodtrack = false;

  //-- |eta| < 2.4
  if(TMath::Abs(eta) > 2.4) isgoodtrack = false;

  if(npixelhit < 3 && TMath::Abs(eta) < 1) isgoodtrack = false;
  if(npixelhit < 2 && TMath::Abs(eta) > 1) isgoodtrack = false;

  if(eta > -1 && eta < 0.4)
    if(phi > 2.4 && phi < 2.8) 
      isgoodtrack = false;

  return(isgoodtrack);
}

//-- good vertex decision 
bool RhoCommon::GetGoodVertex15(bool isfake,bool isvalid,double vz,double rho,int nGoodTrack) {

  bool isgoodvertex = true;

  //-- not fake
  if(isfake) isgoodvertex = false;

  //-- valid
  if(!isvalid) isgoodvertex = false;
  
  //-- rho <= 0.2
  if(rho > 0.2) isgoodvertex = false; 

  //-- |vz| <= 15
  if(TMath::Abs(vz) > 15) isgoodvertex  = false;

  //-- nGoodTrack >= 2
  if(nGoodTrack < 2) isgoodvertex = false;

  return(isgoodvertex);
}

//-- good vertex decision
bool RhoCommon::GetGoodVertex20(bool isfake,bool isvalid,double vz,double rho,int nGoodTrack) {

  bool isgoodvertex = true;

  //-- not fake                                                                                                         
  if(isfake) isgoodvertex = false;

  //-- valid                                                                                                         
  if(!isvalid) isgoodvertex = false;

  //-- rho <= 0.2                                                                                                  
  if(rho > 0.2) isgoodvertex = false;

  //-- |vz| <= 20                                                                                
  if(TMath::Abs(vz) > 20) isgoodvertex= false;

  //-- nGoodTrack >= 2                                                                                    
  if(nGoodTrack < 2) isgoodvertex = false;

  return(isgoodvertex);
}

//-- good particle decision
bool RhoCommon::GetGoodParticle(int status, int charge, double eta, double pT) {

  bool isgoodparticle = true;

  //-- status one
  if(status != 1) isgoodparticle = false;

  //-- charged particle
  if(charge == 0) isgoodparticle = false;

  //-- |eta| <= 2.4
  if(TMath::Abs(eta) > 2.4) isgoodparticle = false;

  //-- pT > 0.5 
  if(pT < 0.5)  isgoodparticle = false;

  return(isgoodparticle);
}

//-- compute correlation coefficient                                                                                    
void RhoCommon::GetCorrelation(TH1D* hnFwd[n_delta_eta], TH1D* hnBwd[n_delta_eta], TH1D* hnFwdBwd[n_delta_eta], TH1D* hrho) {

  cout<<endl;
  //-- cout<<"histo "<<hrho->GetName()<<" with title "<<hrho->GetTitle()<<" will be filled"<<endl<<endl;

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

    hrho->SetBinContent(i+1,rho[i]);
    hrho->SetBinError(i+1,rho_error[i]);
  }                            

}

//-- compute correlation coefficient                                                                                    
void RhoCommon::GetCorrelation(TH2D* hcorrelation[n_delta_eta], TH1D* hrho) {

  cout<<endl;
  //-- cout<<"histo "<<hrho->GetName()<<" with title "<<hrho->GetTitle()<<" will be filled"<<endl<<endl;

  //-- correlation coefficient  
  for(int i = 0; i < n_delta_eta; ++i) {    

    double n = hcorrelation[i]->GetEntries();
    double rho = hcorrelation[i]->GetCorrelationFactor();
    double rho_error= (1-TMath::Power(rho,2))/TMath::Sqrt(n-1);

    double delta_eta = eta_plus_min[i]-eta_minus_max[i];
    cout<<"delta eta = "<<delta_eta<<": rho = "<<rho<<" +/- "<<rho_error<<endl;

    hrho->SetBinContent(i+1,rho);
    hrho->SetBinError(i+1,rho_error);
  }                            

}

//-- compute correlation coefficient                                                                                                                                     
void RhoCommon::GetCorrelation(TH2D* hcorrelation[nHFsel][n_delta_eta], TH1D* hrho[nHFsel]) {

  //-- correlation coefficient                                                                                                                                            
  for(int isel = 0; isel < nHFsel; isel++) {											   
    
    //-- cout<<"histo "<<hrho[isel]->GetName()<<" - "<<hrho[isel]->GetTitle()<<":"<<endl<<endl;

    for(int ieta = 0; ieta < n_delta_eta; ++ieta) {										   
      
      double n = hcorrelation[isel][ieta]->GetEntries();
      double rho = hcorrelation[isel][ieta]->GetCorrelationFactor();
      double rho_error= (1-TMath::Power(rho,2))/TMath::Sqrt(n-1);
      
      //-- double delta_eta = eta_plus_min[ieta]-eta_minus_max[ieta];
      //-- cout<<"delta eta = "<<delta_eta<<": rho = "<<rho<<" +/- "<<rho_error<<endl;
      
      hrho[isel]->SetBinContent(ieta+1,rho);
      hrho[isel]->SetBinError(ieta+1,rho_error);
    }    
    //-- cout<<endl;
  }

}

//-- compute histo ratio
void RhoCommon::ComputeHistoRatio(TH1D* hup[nHFsel], TH1D* hdown[nHFsel], TH1D* hratio[nHFsel]) {

  for(int isel = 0; isel < nHFsel; isel++) {

    for(int ibin = 1; ibin <= hup[isel]->GetNbinsX(); ibin++) {

      double up = hup[isel]->GetBinContent(ibin);
      double down = hdown[isel]->GetBinContent(ibin);

      double error_up = hup[isel]->GetBinError(ibin);
      double error_down = hdown[isel]->GetBinError(ibin);

      double ratio = 0;
      double error = 0;
      Common::ComputeRatio(up,down,error_up,error_down,ratio,error);

      hratio[isel]->SetBinContent(ibin,ratio);
      hratio[isel]->SetBinError(ibin,error);
    }
  }

}



//-- retrieve branches
void RhoCommon::GetBranches(TTree* tree,bool isData) {

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
}

//-- set branch address
void RhoCommon::SetBranchAddresses(TTree* tree,bool isData) {

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
}

//-- get tree entries
void RhoCommon::GetTreeEntries(int ievt, bool isData) {

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
}

//-- create histo
void RhoCommon::CreateHistoSelection() {

  //-- selection
  hselection_data = MakeHisto("hselection_data","Event Selection Data","step","N evts",10,0.5,10.5);
  hselection_moca_reco = MakeHisto("hselection_moca_reco","Event Selection MC reco","step","N evts",7,0.5,7.5);
  hselection_moca_gen = MakeHisto("hselection_moca_gen","Event Selection MC gen","step","N evts",6,0.5,6.5);  
  cout<<endl<<endl;
}

void RhoCommon::CreateHistoHFTower() {
  
  //-- HF tower
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    label = "hnHFTower_" + TString::Format("%d",ihf+1);
    title = "HF tower multiplicity for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hnHFTower[ihf] = Common::MakeHisto(label,title,"N HF towers","N evts",500,0,500);                         

    label = "hHFTowerEtot_" + TString::Format("%d",ihf+1);						   
    title = "HF tower total energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFTowerEtot[ihf] = Common::MakeHisto(label,title,"HF tower Etot","N evts",500,0,5000);
    
    label = "hHFTowerLeadingEnergy_" + TString::Format("%d",ihf+1);				     
    title = "HF tower leading energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFTowerLeadingEnergy[ihf] = Common::MakeHisto(label,title,"HF tower leading energy","N evts",500,0,500);

    label = "hHFTowerEnergy_" + TString::Format("%d",ihf+1);						   
    title = "HF tower energy for E tower > " + TString::Format("%3.1f",EHFcut[ihf]) + " GeV";
    hHFTowerEnergy[ihf] = Common::MakeHisto(label,title,"HF tower energy","N HF towers",500,0,500);

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
  
}

void RhoCommon::CreateHistoHFplusTower() {

  //-- HF plus tower     
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

}

void RhoCommon::CreateHistoHFminusTower() { 

  //-- HF minus tower
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
}

void RhoCommon::CreateHistoTrack() {

  //-- track distributions exactly one good vertex
  hTrackNum = MakeHisto("hTrackNum","track multiplicity","N tracks","N evts",500,0,500);
  hTrackPt = MakeHisto("hTrackPt","track pT","track pT","N tracks",150,0,15);
  hTrackEta = MakeHisto("hTrackEta","track eta","track #eta","N tracks",120,-3,3);

  hTrackPhi = MakeHisto("hTrackPhi","track phi","track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());
  hTrackPhiCentral = MakeHisto("hTrackPhiCentral","track phi |eta| < 1","track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());
  hTrackPhiFwd = MakeHisto("hTrackPhiFwd","track phi |eta| > 1","track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());

  hTrackEtaPhi = MakeHisto("hTrackEtaPhi","track eta - phi","track #eta","track #phi",240,-3,3,200,-TMath::Pi(),TMath::Pi());

  hTrackQuality = MakeHisto("hTrackQuality","track quality","track quality","N tracks",3,-0.5,2.5);

  for(int ieta = 0; ieta < neta_edge-1; ieta++) {
    label = "hTrackPtEta_" + TString::Format("%d",ieta+1);
    title = " track pT for " + TString::Format("%3.1f",eta_edge[ieta]) + " < |eta| < " + TString::Format("%3.1f",eta_edge[ieta+1]);
    hTrackPtEta[ieta] = MakeHisto(label,title,"track pT","N tracks",150,0,15);
  }

  hTrackPtError = MakeHisto("hTrackPtError","track pT error/pT","track pT error/pT","N tracks",100,0,10);
  hTrackEtaError = MakeHisto("hTrackEtaError","track eta error","track #eta error", "N tracks", 100,0,0.02);
  hTrackPhiError = MakeHisto("hTrackPhiError","track phi error","track #phi error","N tracks",100,0,0.02);
  
  hTrackPtErrorPt = MakeProfile("hTrackPtErrorPt","track pT error/pT versus pT","track pT","track pT error/pT",100,0,5);
  hTrackPtErrorEta = MakeProfile("hTrackPtErrorEta","track pT error/pT versus #eta","track #eta","track pT error/pT",100,-2.5,2.5);
  hTrackPtErrorNhit = MakeProfile("hTrackPtErrorNhit","track pT error/pT versus Nhit","track Nhit","track pT error/pT",51,-0.5,50.5);

  hTrackdxy = MakeHisto("hTrackdxy","track dxy","track dxy/sigmaxy","N tracks",100,0,10);
  hTrackdz = MakeHisto("hTrackdz","track dz","track dz/sigmaz","N tracks",100,0,10);

  hTrackdxyEta = MakeProfile("hTrackdxyEta","track dxy/sigmaxy versus eta","track eta","track dxy/sigmaxy",120,-3,3);
  hTrackdzEta = MakeProfile("hTrackdzEta","track dz/sigmaz versus eta","track eta","track dz/sigmaz",120,-3,3);       

  hTrackdxyNhit = MakeProfile("hTrackdxyNhit","track dxy/sigmaxy versus Nhit","track Nhit","track dxy/sigmaxy",51,-0.5,50.5);
  hTrackdzNhit = MakeProfile("hTrackdzNhit","track dz/sigmaz versus Nhit","track Nhit","track dz/sigmaz",51,-0.5,50.5);       
  
  hTrackdxyEta2D = MakeHisto("hTrackdxyEta2D","track dxy/sigmaxy versus eta","track eta","track dxy/sigmaxy",240,-3,3,200,-10,10);
  hTrackdzEta2D = MakeHisto("hTrackdzEta2D","track dz/sigmaz versus eta","track eta","track dz/sigmaz",240,-3,3,200,-10,10);

  hTrackdxyNum2D = MakeHisto("hTrackdxyNum2D","track dxy/sigmaxy versus N tracks","N tracks","track dxy/sigmaxy",200,0,100,200,-10,10);
  hTrackdzNum2D = MakeHisto("hTrackdzNum2D","track dz/sigmaz versus N tracks","N tracks","track dz/sigmaz",200,0,100,200,-10,10);

  hTrackchi2 = MakeHisto("hTrackchi2","track chi2","track chi2","N tracks",100,0,100);
  hTrackndof = MakeHisto("hTrackndof","track ndof","track ndof","N tracks",100,0,100);
  hTrackchi2ndof = MakeHisto("hTrackchi2ndof","track chi2ndof","track chi2ndof","N tracks",100,0,10);

  hTrackValidHit = MakeHisto("hTrackValidHit","track ValidHit","track ValidHit","N tracks",101,-0.5,100.5);
  hTrackValidPixelHit = MakeHisto("hTrackValidPixelHit","track ValidPixelHit","track ValidPixelHit","N tracks",101,-0.5,100.5);
  hTrackValidStripHit = MakeHisto("hTrackValidStripHit","track ValidStripHit","track ValidStripHit","N tracks",101,-0.5,100.5);
 
  hTrackWeight = MakeHisto("hTrackWeight","track weight","track weight","N tracks",100,0,1);

  hTrackValidHitEta = MakeProfile("hTrackValidHitEta","track ValidHit versus eta","eta","track ValidHit",120,-3,3);
  hTrackValidPixelHitEta = MakeProfile("hTrackValidPixelHitEta","track ValidPixelHit versus eta","eta","track ValidPixelHit",120,-3,3);
  hTrackValidStripHitEta = MakeProfile("hTrackValidStripHitEta","track ValidStripHit versus eta","eta","track ValidStripHit",120,-3,3);

  hLeadingTrackPt = MakeHisto("hLeadingTrackPt","leading track pT","leading track pT","N tracks",150,0,15);
  hLeadingTrackEta = MakeHisto("hLeadingTrackEta","leading track eta","leading track #eta","N tracks",120,-3,3);
  hLeadingTrackPhi = MakeHisto("hLeadingTrackPhi","leading track phi","leading track #phi","N tracks",100,-TMath::Pi(),TMath::Pi());

  cout<<endl<<endl;
}

void RhoCommon::CreateHistoVertex() {

  //-- vertex distributions exactly one good vertex
  hVx = MakeHisto("hVx","vertex x ","vertex x [cm]","N evts",800,-0.2,0.2);
  hVy = MakeHisto("hVy","vertex y ","vertex y [cm]","N evts",800,-0.2,0.2);
  hVz = MakeHisto("hVz","vertex z","vertex z [cm]","N evts",210,-21.5,20.5);
  hVxy = MakeHisto("hVxy","vertex y versus x","vertex x [cm]","vertex y [cm]",200,-0.2,0.2,200,-0.2,0.2);
  hVrho = MakeHisto("hVrho","vertex rho ","vertex #rho [cm]","N evts",10000,0,1);
  
  hVchi2 = MakeHisto("hVchi2","vertex chi2","vertex #chi^{2}","N evts",500,0,100);
  hVndof = MakeHisto("hVndof","vertex ndof","vertex ndof","N evts",100,0,100);
  hVchi2ndof = MakeHisto("hVchi2ndof","vertex chi2/ndof","vertex #chi^{2}/ndof","N evts",200,0,10); 
  
  hVzError = MakeHisto("hVzError","vertex z error","vertex z error [cm]","N evts",200,0,0.2);
  hVzErrorVz = MakeHisto("hVzErrorVz","vertex z error versus vertex z","vertex z [cm]","vertex z error [cm]",100,-15,15,100,0,0.2);
  hVzErrorTrackNum = MakeHisto("hVzErrorTrackNum","vertex z error versus track multiplicity","N tracks","vertex z error [cm]",100,0,100,100,0,0.2);

  hVNum = MakeHisto("hVNum","vertex multiplicity","vertex multiplicity","N events",10,0.5,10.5);
  hVNumGood = MakeHisto("hVNumGood","good vertex multiplicity","good vertex multiplicity","N events",10,0.5,10.5);

  cout<<endl<<endl;
}

void RhoCommon::CreateHistoBS() {

  //-- beam spot distributions
  hBSx = MakeHisto("hBSx","BS x","BS x [cm]","N evts",1000,-0.21,0.19);
  hBSy = MakeHisto("hBSy","BS y","BS y [cm]","N evts",1000,-0.21,0.19);
  hBSz = MakeHisto("hBSz","BS z","BS z [cm]","N evts",1000,-20.05,19.95);
  
  cout<<endl<<endl;
}

void RhoCommon::CreateHistoPU() {

  //-- pileup information
  hnbxint = MakeHisto("hnbxint","in time bunch crossing multiplicity","in time bunch crossing multiplicity","number of events",11,-0.5,10.5);
  hnbxearly = MakeHisto("hnbxearly","early bunch crossing multiplicity","early bunch crossing multiplicity","number of events",11,-0.5,10.5);
  hnbxlate = MakeHisto("hnbxlate","late bunch crossing multiplicity","late bunch crossing multiplicity","number of events",11,-0.5,10.5);
  hnbxtot = MakeHisto("hnbxtot","total bunch crossing multiplicity","total bunch crossing multiplicity","number of events",11,-0.5,10.5);  

  hPUint = MakeHisto("hPUint","1 + in time PU multiplicity","1 + in time PU multiplicity","number of events",11,-0.5,10.5);
  hPUearly = MakeHisto("hPUearly","1 + early PU multiplicity","1 + early PU multiplicity","number of events",11,-0.5,10.5);
  hPUlate = MakeHisto("hPUlate","1 + late PU multiplicity","1 + late PU multiplicity","number of events",11,-0.5,10.5);
  hPUtot = MakeHisto("hPUtot","1 + total PU multiplicity","1 + total PU multiplicity","number of events",11,-0.5,10.5);

  hNumInt = MakeHisto("hNumInt","interactions multiplicity","interactions multiplicity","number of events",100,0,10); 

  hPUintVNum = MakeHisto("hPUintVNum","vertex multi versus 1 + in time PU multi","1 + in time PU multi","vertex multi",11,-0.5,10.5,11,-0.5,10.5);
  hPUearlyVNum = MakeHisto("hPUearlyVNum","vertex multi versus 1 + early PU multi","1 + early PU multi","vertex multi",11,-0.5,10.5,11,-0.5,10.5);
  hPUlateVNum = MakeHisto("hPUlateVNum","vertex multi versus 1 + late PU multi","1 + late PU multi","vertex multi",11,-0.5,10.5,11,-0.5,10.5);
  hPUtotVNum = MakeHisto("hPUtotVNum","vertex multi versus 1 + total PU multi","1 + total PU multi","vertex multi",11,-0.5,10.5,11,-0.5,10.5);      
  hPUintVNum0 = MakeHisto("hPUintVNum0","vertex multi versus 1 + in time PU multi","1 + in time PU multi","vertex multi",11,-0.5,10.5,11,-0.5,10.5);

  cout<<endl<<endl;
}

void RhoCommon::CreateHistoRho() {

  //-- correlation coefficient  							    
  for(int i = 0; i < n_delta_eta; ++i) {						    
    double delta_eta = i*delta_eta_0;							    
    											    
    eta_plus_min[i] = 0.5*delta_eta;							    
    eta_plus_max[i] = 0.5*delta_eta + delta_eta_0;					    
    cout<<"forward interval: ["<<eta_plus_min[i]<<","<<eta_plus_max[i]<<"]"<<endl;     
    											    
    eta_minus_min[i] = -0.5*delta_eta - delta_eta_0;					    
    eta_minus_max[i] = -0.5*delta_eta;							    
    cout<<"backward interval: ["<<eta_minus_min[i]<<","<<eta_minus_max[i]<<"]"<<endl;  
  }
                                                    				    
  cout<<endl;                                                                             

  //-- correlation coefficient  											     
  for(int isel = 0; isel < nHFsel; isel++) {
    for(int ieta = 0; ieta < n_delta_eta; ++ieta) {										   

      label = "hcorrelation_sel_" + TString::Format("%d",isel) + "_eta_" + TString::Format("%d",ieta);							  
      title_fwd = TString::Format("%3.2f",eta_plus_min[ieta]) + " < eta fwd < " + TString::Format("%3.2f",eta_plus_max[ieta]);	   
      title_bwd = TString::Format("%3.2f",eta_minus_min[ieta]) + " < eta bwd < " + TString::Format("%3.2f",eta_minus_max[ieta]);   
      title = HFsel[isel] + " - n bwd versus n fwd - " + title_fwd + " - " + title_bwd;						   
      hcorrelation[isel][ieta] = MakeHisto(label,title,"n fwd","n bwd",40,0,40,40,0,40);                                              

      label = "hcorrelationDown_sel_" + TString::Format("%d",isel) + "_eta_" + TString::Format("%d",ieta);				      
      title_fwd = TString::Format("%3.2f",eta_plus_min[ieta]) + " < eta fwd < " + TString::Format("%3.2f",eta_plus_max[ieta]);	      
      title_bwd = TString::Format("%3.2f",eta_minus_min[ieta]) + " < eta bwd < " + TString::Format("%3.2f",eta_minus_max[ieta]);      
      title = HFsel[isel] + " - n bwd versus n fwd - " + title_fwd + " - " + title_bwd;						      
      hcorrelationDown[isel][ieta] = MakeHisto(label,title,"n fwd","n bwd",40,0,40,40,0,40);                                              

      label = "hcorrelationUp_sel_" + TString::Format("%d",isel) + "_eta_" + TString::Format("%d",ieta);				  
      title_fwd = TString::Format("%3.2f",eta_plus_min[ieta]) + " < eta fwd < " + TString::Format("%3.2f",eta_plus_max[ieta]);	      	  
      title_bwd = TString::Format("%3.2f",eta_minus_min[ieta]) + " < eta bwd < " + TString::Format("%3.2f",eta_minus_max[ieta]);      	  
      title = HFsel[isel] + " - n bwd versus n fwd - " + title_fwd + " - " + title_bwd;						      	  
      hcorrelationUp[isel][ieta] = MakeHisto(label,title,"n fwd","n bwd",40,0,40,40,0,40);                                              

      label = "hcorrelationPU_sel_" + TString::Format("%d",isel) + "_eta_" + TString::Format("%d",ieta);
      title_fwd = TString::Format("%3.2f",eta_plus_min[ieta]) + " < eta fwd < " + TString::Format("%3.2f",eta_plus_max[ieta]);	  
      title_bwd = TString::Format("%3.2f",eta_minus_min[ieta]) + " < eta bwd < " + TString::Format("%3.2f",eta_minus_max[ieta]);  
      title = HFsel[isel] + " - n bwd versus n fwd - " + title_fwd + " - " + title_bwd;
      hcorrelationPU[isel][ieta] = MakeHisto(label,title,"n fwd","n bwd",40,0,40,40,0,40);                                               

      label = "hcorrelationEffi_sel_" + TString::Format("%d",isel) + "_eta_" + TString::Format("%d",ieta);				 
      title_fwd = TString::Format("%3.2f",eta_plus_min[ieta]) + " < eta fwd < " + TString::Format("%3.2f",eta_plus_max[ieta]);	  	 
      title_bwd = TString::Format("%3.2f",eta_minus_min[ieta]) + " < eta bwd < " + TString::Format("%3.2f",eta_minus_max[ieta]);  	 
      title = HFsel[isel] + " - n bwd versus n fwd - " + title_fwd + " - " + title_bwd;							 
      hcorrelationEffi[isel][ieta] = MakeHisto(label,title,"n fwd","n bwd",40,0,40,40,0,40);                                               

      label = "hcorrelationGen_sel_" + TString::Format("%d",isel) + "_eta_" + TString::Format("%d",ieta);				   
      title_fwd = TString::Format("%3.2f",eta_plus_min[ieta]) + " < eta fwd < " + TString::Format("%3.2f",eta_plus_max[ieta]);	  	   
      title_bwd = TString::Format("%3.2f",eta_minus_min[ieta]) + " < eta bwd < " + TString::Format("%3.2f",eta_minus_max[ieta]);  	   
      title = HFsel[isel] + " - n bwd versus n fwd - " + title_fwd + " - " + title_bwd;							   
      hcorrelationGen[isel][ieta] = MakeHisto(label,title,"n fwd","n bwd",40,0,40,40,0,40);                                               

      label = "hcorrelationTheo_sel_" + TString::Format("%d",isel) + "_eta_" + TString::Format("%d",ieta);				  
      title_fwd = TString::Format("%3.2f",eta_plus_min[ieta]) + " < eta fwd < " + TString::Format("%3.2f",eta_plus_max[ieta]);	  	  
      title_bwd = TString::Format("%3.2f",eta_minus_min[ieta]) + " < eta bwd < " + TString::Format("%3.2f",eta_minus_max[ieta]);  	  
      title = HFsel[isel] + " - n bwd versus n fwd - " + title_fwd + " - " + title_bwd;							  
      hcorrelationTheo[isel][ieta] = MakeHisto(label,title,"n fwd","n bwd",40,0,40,40,0,40);                                               
    }
  }

  cout<<endl;

  double betamin = - 0.5*delta_eta_0;
  double betamax = betamin + n_delta_eta*delta_eta_0;

  for(int isel = 0; isel < nHFsel; isel++) {	   					
    
    label = "hrho_sel_" + TString::Format("%d",isel);
    title = HFsel[isel] + " - correlation coefficient";
    hrho[isel] = MakeHisto(label,title,"delta eta","rho fwd-bwd",n_delta_eta,betamin,betamax);
    
    label = "hrhoDown_sel_" + TString::Format("%d",isel);					      
    title = HFsel[isel] + " - correlation coefficient";					      
    hrhoDown[isel] = MakeHisto(label,title,"delta eta","rho fwd-bwd",n_delta_eta,betamin,betamax);

    label = "hrhoUp_sel_" + TString::Format("%d",isel);					      
    title = HFsel[isel] + " - correlation coefficient";					      
    hrhoUp[isel] = MakeHisto(label,title,"delta eta","rho fwd-bwd",n_delta_eta,betamin,betamax);

    label = "hrhoPU_sel_" + TString::Format("%d",isel);								      
    title = HFsel[isel] + " - correlation coefficient";				      
    hrhoPU[isel] = MakeHisto(label,title,"delta eta","rho fwd-bwd",n_delta_eta,betamin,betamax);

    label = "hrhoEffi_sel_" + TString::Format("%d",isel);						
    title = HFsel[isel] + " - correlation coefficient";				      		
    hrhoEffi[isel] = MakeHisto(label,title,"delta eta","rho fwd-bwd",n_delta_eta,betamin,betamax);

    label = "hrhoGen_sel_" + TString::Format("%d",isel);					  
    title = HFsel[isel] + " - correlation coefficient";				      		  
    hrhoGen[isel] = MakeHisto(label,title,"delta eta","rho fwd-bwd",n_delta_eta,betamin,betamax);

    label = "hrhoTheo_sel_" + TString::Format("%d",isel);					 
    title = HFsel[isel] + " - correlation coefficient";				      		 
    hrhoTheo[isel] = MakeHisto(label,title,"delta eta","rho fwd-bwd",n_delta_eta,betamin,betamax);

    label = "hrhoCF_sel_" + TString::Format("%d",isel);
    title = HFsel[isel] + " - correlation coefficient CF";
    hrhoCF[isel] = MakeHisto(label,title,"delta eta","correction factor",n_delta_eta,betamin,betamax);

    label = "hrhoCFPU_sel_" + TString::Format("%d",isel);									      
    title = HFsel[isel] + " - correlation coefficient CF PU";							      
    hrhoCFPU[isel] = MakeHisto(label,title,"delta eta","correction factor",n_delta_eta,betamin,betamax);
  }

  cout<<endl<<endl;
}

//-- check histo
void RhoCommon::CheckHistoSelection() {

  //-- selection
  CheckHisto(hselection_data);
  CheckHisto(hselection_moca_reco);
  CheckHisto(hselection_moca_gen);
  cout<<endl<<endl;
}

void RhoCommon::CheckHistoHFTower() {

  //-- HF tower   
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    CheckHisto(hnHFTower[ihf]);
    CheckHisto(hHFTowerEtot[ihf]);
    CheckHisto(hHFTowerLeadingEnergy[ihf]);

    CheckHisto(hHFTowerEnergy[ihf]);
    CheckHisto(hHFTowerEta[ihf]);
    CheckHisto(hHFTowerPhi[ihf]);

    CheckHisto(hHFTowerEtaPhi[ihf]);
  }

  cout<<endl<<endl;
}
  
void RhoCommon::CheckHistoHFplusTower() {

  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    CheckHisto(hnHFplusTower[ihf]);   
    CheckHisto(hHFplusTowerEtot[ihf]);   
    CheckHisto(hHFplusTowerLeadingEnergy[ihf]);

    CheckHisto(hHFplusTowerEnergy[ihf]);   
    CheckHisto(hHFplusTowerEta[ihf]);   
    CheckHisto(hHFplusTowerPhi[ihf]);    

    CheckHisto(hHFplusTowerEtaPhi[ihf]);
  }
  
  cout<<endl<<endl;
}

void RhoCommon::CheckHistoHFminusTower() {

  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    CheckHisto(hnHFminusTower[ihf]);       
    CheckHisto(hHFminusTowerEtot[ihf]);       
    CheckHisto(hHFminusTowerLeadingEnergy[ihf]);

    CheckHisto(hHFminusTowerEnergy[ihf]);       
    CheckHisto(hHFminusTowerEta[ihf]);       
    CheckHisto(hHFminusTowerPhi[ihf]);    

    CheckHisto(hHFminusTowerEtaPhi[ihf]);
 }
  
  cout<<endl<<endl;
}

void RhoCommon::CheckHistoTrack() {

  //-- track distributions exactly one good vertex
  CheckHisto(hTrackNum);
  CheckHisto(hTrackQuality);
  CheckHisto(hTrackPt);
  CheckHisto(hTrackEta);
  
  CheckHisto(hTrackPhi);
  CheckHisto(hTrackPhiCentral);
  CheckHisto(hTrackPhiFwd);
  
  CheckHisto(hTrackEtaPhi);
  
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

  cout<<endl<<endl;
}

void RhoCommon::CheckHistoVertex() {
  
  //-- vertex distributions exactly one good vertex
  CheckHisto(hVx);
  CheckHisto(hVy);
  CheckHisto(hVxy);
  CheckHisto(hVz);
  CheckHisto(hVrho);
  CheckHisto(hVchi2);
  CheckHisto(hVndof);
  CheckHisto(hVchi2ndof);
  CheckHisto(hVzError);
  CheckHisto(hVzErrorVz);
  CheckHisto(hVzErrorTrackNum);
  CheckHisto(hVNum);
  CheckHisto(hVNumGood);

  cout<<endl<<endl;
}

void RhoCommon::CheckHistoBS() {

  //-- beam spot
  CheckHisto(hBSx);
  CheckHisto(hBSy);
  CheckHisto(hBSz);

  cout<<endl<<endl;
}

void RhoCommon::CheckHistoPU(bool isData) {

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

    cout<<endl<<endl;
  }
              
}           

void RhoCommon::CheckHistoRho(bool isData) {

  for(int isel = 0; isel < nHFsel; isel++) {	   						   
    for(int ieta = 0; ieta < n_delta_eta; ++ieta) {
      CheckHisto(hcorrelation[isel][ieta]);

      CheckHisto(hcorrelationDown[isel][ieta]);
      CheckHisto(hcorrelationUp[isel][ieta]);

      CheckHisto(hcorrelationPU[isel][ieta]);
      CheckHisto(hcorrelationEffi[isel][ieta]);

      if(isData) continue;	    
				    
      CheckHisto(hcorrelationGen[isel][ieta]);   	    
      CheckHisto(hcorrelationTheo[isel][ieta]);				
    }
  }

  cout<<endl;

  for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrho[isel]);

  for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoDown[isel]);
  for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoUp[isel]);

  for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoPU[isel]);
  for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoEffi[isel]);

  if(!isData) for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoGen[isel]);
  if(!isData) for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoTheo[isel]);

  if(!isData) for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoCF[isel]);
  if(!isData) for(int isel = 0; isel < nHFsel; isel++) CheckHisto(hrhoCFPU[isel]);

  cout<<endl<<endl;   
}

//-- write histo
void RhoCommon::WriteHistoSelection() {

  //-- selection   
  hselection_data->Write();   
  hselection_moca_reco->Write();
  hselection_moca_gen->Write();
}

void RhoCommon::WriteHistoHFTower() {

  //-- HF tower      
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    hnHFTower[ihf]->Write();
    hHFTowerEtot[ihf]->Write();
    hHFTowerLeadingEnergy[ihf]->Write();

    hHFTowerEnergy[ihf]->Write();
    hHFTowerEta[ihf]->Write();
    hHFTowerPhi[ihf]->Write();    

    hHFTowerEtaPhi[ihf]->Write();
  }
}

void RhoCommon::WriteHistoHFplusTower() {

  //-- HF plus tower      
  for(int ihf = 0; ihf < nHFcut; ++ihf) {    
    hnHFplusTower[ihf]->Write();
    hHFplusTowerEtot[ihf]->Write();   
    hHFplusTowerLeadingEnergy[ihf]->Write();

    hHFplusTowerEnergy[ihf]->Write();
    hHFplusTowerEta[ihf]->Write();   
    hHFplusTowerPhi[ihf]->Write();   

    hHFplusTowerEtaPhi[ihf]->Write();
  }
}

void RhoCommon::WriteHistoHFminusTower() {

  //-- HF minus tower      
  for(int ihf = 0; ihf < nHFcut; ++ihf) {
    hnHFminusTower[ihf]->Write();       
    hHFminusTowerEtot[ihf]->Write();  
    hHFminusTowerLeadingEnergy[ihf]->Write();

    hHFminusTowerEnergy[ihf]->Write();       
    hHFminusTowerEta[ihf]->Write();       
    hHFminusTowerPhi[ihf]->Write();   
    
    hHFminusTowerEtaPhi[ihf]->Write();
  }
}

void RhoCommon::WriteHistoTrack() {
  
  //-- track distributions exactly one good vertex
  hTrackNum->Write(); 
  hTrackQuality->Write();
  hTrackPt->Write();
  hTrackEta->Write();
  
  hTrackPhi->Write();
  hTrackPhiCentral->Write();
  hTrackPhiFwd->Write();
  
  hTrackEtaPhi->Write();
  
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
}

void RhoCommon::WriteHistoVertex() {

  //-- vertex distributions exactly one good vertex
  hVx->Write();
  hVy->Write();
  hVxy->Write();
  hVz->Write();
  hVrho->Write();

  hVchi2->Write();
  hVndof->Write();
  hVchi2ndof->Write();

  hVzError->Write();
  hVzErrorVz->Write();
  hVzErrorTrackNum->Write();

  hVNum->Write();
  hVNumGood->Write();
}

void RhoCommon::WriteHistoBS() {

  //-- beam spot   
  hBSx->Write();
  hBSy->Write();
  hBSz->Write();
}

void RhoCommon::WriteHistoPU(bool isData) {

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
}

void RhoCommon::WriteHistoRho(bool isData) {

  //-- correlation coefficient    
  for(int isel = 0; isel < nHFsel; isel++) {	      
    for(int ieta = 0; ieta < n_delta_eta; ++ieta) {   
						      
      hcorrelation[isel][ieta]->Write();

      hcorrelationDown[isel][ieta]->Write();
      hcorrelationUp[isel][ieta]->Write();

      hcorrelationPU[isel][ieta]->Write();						      
      hcorrelationEffi[isel][ieta]->Write();						      

      if(isData) continue;	 
				 
      hcorrelationGen[isel][ieta]->Write();        
      hcorrelationTheo[isel][ieta]->Write();			
    }
  }

  for(int isel = 0; isel < nHFsel; isel++) hrho[isel]->Write();

  for(int isel = 0; isel < nHFsel; isel++) hrhoDown[isel]->Write();
  for(int isel = 0; isel < nHFsel; isel++) hrhoUp[isel]->Write();

  for(int isel = 0; isel < nHFsel; isel++) hrhoPU[isel]->Write();
  for(int isel = 0; isel < nHFsel; isel++) hrhoEffi[isel]->Write();

  if(!isData) for(int isel = 0; isel < nHFsel; isel++) hrhoGen[isel]->Write();
  if(!isData) for(int isel = 0; isel < nHFsel; isel++) hrhoTheo[isel]->Write();

  if(!isData) for(int isel = 0; isel < nHFsel; isel++) hrhoCF[isel]->Write();
  if(!isData) for(int isel = 0; isel < nHFsel; isel++) hrhoCFPU[isel]->Write();
}                                        

