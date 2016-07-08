#include "RhoAnalyzer.h"

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

#define CheckHFsel 0
#define CheckPUWeight 0
#define CheckGoodTrackVertex 0
#define CheckVertexWeight 0
#define CheckEffi 0
#define CheckScalarPtSum 0

RhoAnalyzer::RhoAnalyzer() { }

RhoAnalyzer::~RhoAnalyzer() { }

void RhoAnalyzer::Loop(TString inputdir,TObjArray* filelist,TString type, TString selection, TString file_name) {

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

  TString MCname[3];
  MCname[0] = "PYTHIA8";
  MCname[1] = "HERWIGPP";
  MCname[2] = "EPOS";

  bool AllMC[3];
  for(int i = 0; i < 3; i++) AllMC[i] = false;
  if(isPYTHIA8 == true) AllMC[0] = true;
  if(isHERWIGPP == true) AllMC[1] = true;
  if(isEPOS == true) AllMC[2] = true;

  bool pu_reweight = false;
  bool pixel_cut = false;
  bool vz_reweight = false; 

  double HFscaling = 0;
				 
  if(strcmp(selection,"no_pu_reweight") == 0) {  
   pu_reweight = false; 		       
   pixel_cut = true;   		       
   vz_reweight = true; 		       

   if(isData) HFscaling = 1;
   if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }                                            

  if(strcmp(selection,"no_pixel_cut") == 0) {
    pu_reweight = true; 
    pixel_cut = false;   
    vz_reweight = true; 

    if(isData) HFscaling = 1;				
    if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }

  if(strcmp(selection,"no_vz_reweight") == 0) {
   pu_reweight = true; 
   pixel_cut = true;   
   vz_reweight = false; 

   if(isData) HFscaling = 1;				   
   if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }

  if(strcmp(selection,"nominal") == 0) {
    pu_reweight = true;
    pixel_cut = true;
    vz_reweight = true;

   if(isData) HFscaling = 1;				
   if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }                                               

  if(strcmp(selection,"HFup") == 0) {
    pu_reweight = true;
    pixel_cut = true;
    vz_reweight = true;

    if(isData) HFscaling = 1.20;				
    if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }
 
  if(strcmp(selection,"HFdown")  == 0) {
    pu_reweight = true;
    pixel_cut = true;
    vz_reweight = true;

    if(isData) HFscaling = 0.80;			
    if(!isData) HFscaling = 0.895; //-- 1/1.117	        
  }

  cout<<endl;
  cout<<"type = "<<type<<endl;
  cout<<"selection = "<<selection<<endl<<endl;
  cout<<"pu reweight = "<<pu_reweight<<endl;
  cout<<"pixel cut = "<<pixel_cut<<endl;
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

  long double Nmoca_gen_sel = 0;

  double nSelTrack100 = 0; 
  double nSelTrack96 = 0;  

  double nSelTrack100_tot = 0; 
  double nSelTrack96_tot = 0;  
  
  //-- Create histo
  CreateHistoSelection();

  CreateHistoHFTower();
  CreateHistoHFplusTower();
  CreateHistoHFminusTower();

  CreateHistoTrack();
  CreateHistoVertex();
  CreateHistoBS();  

  CreateHistoPU();                   
  CreateHistoRho();

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
      cout<<"Error in RhoAnalyzer: could not open file "<<temp_itfile->GetString()<<endl;
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

  TRandom3* random = new TRandom3();
  double trackeffi_cut = 0.96;

  while((itfile = (TObjString*)next()) && (file_nb < file_nb_max || file_nb_max == -1)) {
    
    file_nb++;
    
    filename.Clear();
    filename = itfile->GetString();
    
    cout<<endl<<"open file "<<file_nb<<" with name: "<<itfile->GetString()<<endl;
    
    TFile* file = TFile::Open(inputdir+itfile->GetString(),"READ");
    
    if (!file) { 
      cout<<"Error in RhoAnalyzer: could not open file "<<itfile->GetString()<<endl; 
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

      //-- rho with respect to the beam spot
      double rhoBS[nVertex];
      double VxBS[nVertex];
      double VyBS[nVertex];

      for (int ivertex = 0; ivertex < nVertex; ++ivertex) {	
	double deltax = Vx->at(ivertex) - BSx;
	double deltay = Vy->at(ivertex) - BSy;
	VxBS[ivertex] = deltax;
	VyBS[ivertex] = deltay;
	rhoBS[ivertex] = TMath::Sqrt(deltax*deltax+deltay*deltay);
      }
      
      //-- good track flag - initialized for each event 
      std::vector < std::vector<bool> > isGoodTrack;

      //-- good track number
      int nGoodTrack[nVertex];
      for (int ivertex = 0; ivertex < nVertex; ivertex++) nGoodTrack[ivertex] = 0;

      //-- good vertex flag
      bool isGoodVertex15[nVertex];
      for (int ivertex = 0; ivertex < nVertex; ivertex++) isGoodVertex15[ivertex] = false;

      bool isGoodVertex20[nVertex];							 
      for (int ivertex = 0; ivertex < nVertex; ivertex++) isGoodVertex20[ivertex] = false;

      //-- vertex scalar pT sum								  
      double ScalarPtSum[nVertex];
      for (int ivertex = 0; ivertex < nVertex; ivertex++) ScalarPtSum[ivertex] = 0; 

      double ScalarPtSum_max = 0;
      int iGoodVertex15_max = -1;

      //-- good vertex number
      int nGoodVertex15 = 0;
      int nGoodVertex20 = 0;

      //-- loop over vertices
      for (int ivertex = 0; ivertex < nVertex; ++ivertex) {	
	
	//-- initialized for each vertex
	std::vector<bool> isGoodTrackVertex;

	//-- loop over associated tracks
	for (int itrack = 0; itrack < VnTrack->at(ivertex); itrack++) {  
	  
	  int quality = VTrackQuality->at(ivertex).at(itrack);
	  		  
	  double pTerror = 1000;
	  if(VTrackPt->at(ivertex).at(itrack) != 0) 
	    pTerror = 100*VTrackPtError->at(ivertex).at(itrack)/VTrackPt->at(ivertex).at(itrack);
	
	  double Sigxy = 1000;
	  if(VTrackdxyError->at(ivertex).at(itrack) != 0)
	    Sigxy = VTrackdxy->at(ivertex).at(itrack)/VTrackdxyError->at(ivertex).at(itrack);
	
	  double Sigz = 1000;	  
	  if(VTrackdzError->at(ivertex).at(itrack) != 0) 
	    Sigz = VTrackdz->at(ivertex).at(itrack)/VTrackdzError->at(ivertex).at(itrack);

	  int npixelhit = VTrackValidPixelHit->at(ivertex).at(itrack);

	  double eta = VTrackEta->at(ivertex).at(itrack);
	  
	  double phi = VTrackPhi->at(ivertex).at(itrack);

	  //-- good track decision
	  isGoodTrackVertex.push_back(GetGoodTrack(quality,pTerror,Sigxy,Sigz,npixelhit,eta,phi));
	  if(isGoodTrackVertex.at(itrack)) nGoodTrack[ivertex]++;

	  //-- vertex scalar pT sum
	  if(isGoodTrackVertex.at(itrack)) ScalarPtSum[ivertex]+=VTrackPt->at(ivertex).at(itrack);

	} //-- end loop over associated tracks

	isGoodTrack.push_back(isGoodTrackVertex);
	isGoodTrackVertex.clear();

	bool isfake =  VisFake->at(ivertex);
	bool isvalid = VisValid->at(ivertex);

	double vz = Vz->at(ivertex);
	double rho = rhoBS[ivertex];

	isGoodVertex15[ivertex] = GetGoodVertex15(isfake,isvalid,vz,rho,nGoodTrack[ivertex]);
	if(isGoodVertex15[ivertex]) nGoodVertex15++;       
                                       
	if(isGoodVertex15[ivertex] && ScalarPtSum[ivertex] > ScalarPtSum_max) {
	  ScalarPtSum_max = ScalarPtSum[ivertex];
	  iGoodVertex15_max = ivertex;
	}
	
	isGoodVertex20[ivertex] = GetGoodVertex20(isfake,isvalid,vz,rho,nGoodTrack[ivertex]);  
	if(isGoodVertex20[ivertex]) nGoodVertex20++;                                              

      } //-- end loop over vertices

      if(CheckGoodTrackVertex) {
	cout<<"-----------------------"<<endl<<endl;
	cout<<"n vertex from nVertex = "<<nVertex<<endl;
	cout<<"n vertex from isGoodTrack = "<<isGoodTrack.size()<<endl;
	cout<<"n good vertex 20 = "<<nGoodVertex20<<endl;
	cout<<"n good vertex 15 = "<<nGoodVertex15<<endl<<endl;
	
	for (int ivertex = 0; ivertex < nVertex; ++ivertex) {	
	  cout<<"vertex "<<ivertex+1<<endl;
	  cout<<"n track from VnTrack "<<VnTrack->at(ivertex)<<endl;
	  cout<<"n track from isGoodTrack: "<<isGoodTrack.at(ivertex).size()<<endl;
	  cout<<"n good track = "<<nGoodTrack[ivertex]<<endl<<endl;
	}
	getchar();
      }  

      //-- vertex information - at least one good vertex
      bool at_least_one_gv_15 = false;
      if(nGoodVertex15 >= 1) at_least_one_gv_15 = true;

      //-- vertex information - exactly one good vertex
      bool exactly_one_gv_15 = false;					      
      if(nGoodVertex15 == 1 && nVertex == 1) exactly_one_gv_15 = true;

      bool exactly_one_gv_20 = false;
      if(nGoodVertex20 == 1 && nVertex == 1) exactly_one_gv_20 = true;

      //-- vertex weight montecarlo 
      weight_vertex = 1; 
      int Niteration_vertex = 2;
      if(!isData && exactly_one_gv_20 && vz_reweight) weight_vertex = GetVertexWeight(Vz->at(0),AllMC,Niteration_vertex); //-- otherwise set to 1

      if(CheckVertexWeight) {                       
	cout<<"n vertex = "<<nVertex<<" - n good vertex 20 = "<<nGoodVertex20<<endl;
	cout<<"vz = "<<Vz->at(0)<<" cm - vertex weight = "<<weight_vertex<<endl;
      }

      //-- pu weight monte carlo
      weight_pu = 1;
      if(!isData && pu_reweight) weight_pu = weight_pu_array[int(PUint)];

      //-- global weight
      weight = weight_pu*weight_vertex;

      if(CheckPUWeight) {
	cout<<"in time PU multiplicity = "<<PUint<<" - weight PU = "<<weight_pu<<endl;
	cout<<"weight = "<<weight<<endl;
	//getchar();
      }

      //-- vertex reconstruction efficiency
      hPUintVNum0->Fill(1+PUint,nGoodVertex15,weight);

      //-- number of events 
      Nevt_tot+=weight;
      Nevt_tot_no_weight+=1;
      Nevt_tot_no_pu_weight+=weight_vertex;
      Nevt_tot_no_vz_weight+=weight_pu; 

      if(isData) hselection_data->Fill(1,weight);
      if(!isData) hselection_moca_reco->Fill(1,weight);
      if(!isData) hselection_moca_gen->Fill(1,weight);

      if((int (Nevt_tot_no_weight))%20000 == 0) {
	cout<<"total number of events done - no weight = "<<Nevt_tot_no_weight<<" ("<<100*Nevt_tot_no_weight/Nevt_all_files<<"%)"<<endl;
	cout<<"total number of events done - global weight = "<<Nevt_tot<<endl<<endl;
	if(isData) cout<<"number of selected events in data - global weight = "<<Ndata_sel<<endl<<endl;
	if(!isData) cout<<"number of selected events in mc at reco level - global weight = "<<Nmoca_reco_sel<<endl;
	if(!isData) cout<<"number of selected events in mc at gen level - global weight = "<<Nmoca_gen_sel<<endl<<endl;
      }       

      //-- Filter on data       
      bool run_selection = false;
      bool L1TT_selection = false;
      bool HLT_selection = false;
      bool data_sel_at_least_one_gv = false;
      bool data_sel = false;
      
      if(isData) {
	
	//-- filter on run and lumi section 	
	if(Run == 251721 && (LumiSection >= 123 && LumiSection <= 244)) run_selection = true;
	
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

	if(run_selection && L1TT_selection && HLT_selection && at_least_one_gv_15) {
	  hselection_data->Fill(5,weight);
	  data_sel_at_least_one_gv = true;
	}

	if(run_selection && L1TT_selection && HLT_selection && exactly_one_gv_15) {
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
	if(at_least_one_gv_15) {
	  hselection_moca_reco->Fill(2,weight);
	  moca_reco_sel_at_least_one_gv = true;
	}
	
	if(exactly_one_gv_15) {
	  hselection_moca_reco->Fill(3,weight);
	  moca_reco_sel = true;

	  Nmoca_reco_sel+=weight;									  
	  Nmoca_reco_sel_no_weight+=1;
	  Nmoca_reco_sel_no_pu_weight+=weight_vertex;
	  Nmoca_reco_sel_no_vz_weight+=weight_pu;
	}
      } //-- end loop on moca reco

      //-- Filter on moca gen
      bool moca_gen_sel = false;

      vector<bool> isGoodParticle;

      int nGoodParticle = 0;
      int ngenparticle = 0;

      bool activity_gen_plus = false; 		       
      bool activity_gen_minus = false;		       
      							       
      bool HFgensel[4]; //-- inclusive - inelastic - NSD -SD 
      for(int i = 0; i < nHFsel; i++) HFgensel[i] = false;   

      if(!isData) {
      
        ngenparticle = ParticleEnergy->size();

	//-- loop over generated particles        
	for(int igen = 0; igen < ngenparticle; ++igen) {

	  int status = ParticleStatus->at(igen);  
	  int charge = ParticleCharge->at(igen);  
	  
	  double eta = ParticleEta->at(igen);	  
	  double pT = ParticlePt->at(igen);       
	  double energy = ParticleEnergy->at(igen);

	  //-- forward selection
	  if(GetHFActivity(status,energy,eta) && eta > 0) activity_gen_plus = true;
	  if(GetHFActivity(status,energy,eta) && eta < 0) activity_gen_minus = true;

	  //-- central selection
	  isGoodParticle.push_back(GetGoodParticle(status,charge,eta,pT));
	  if(isGoodParticle.at(igen)) nGoodParticle++;
	} //-- end loop over generated particles

	//-- inclusive: no HF						   	   
	HFgensel[0] = true;						   	   
      	
	//-- inelastic: HF OR						   	   
	if(activity_gen_plus || activity_gen_minus) HFgensel[1] = true;	   
      									   	   
	//-- NSD: HF AND						   		   
	if(activity_gen_plus && activity_gen_minus) HFgensel[2] = true;	   
      									   	   
	//-- SD: HF XOR							   	   
	if(activity_gen_plus && !activity_gen_minus) HFgensel[3] = true;	   
	if(!activity_gen_plus && activity_gen_minus) HFgensel[3] = true;	   

	//-- at least one charged particle with pT > 0.5 GeV in |eta| < 2.4
	if(nGoodParticle >= 1) moca_gen_sel = true;
	if(nGoodParticle >= 1) Nmoca_gen_sel+=weight;	
      }

      //-- Apply filter - moca gen
      if(!isData && moca_gen_sel) {
	
	hselection_moca_gen->Fill(2,weight);
	
	double nFwdGen[nHFsel][n_delta_eta];
	double nBwdGen[nHFsel][n_delta_eta];

	for(int isel = 0; isel < nHFsel; isel++) {
	  for(int ieta = 0; ieta < n_delta_eta; ++ieta) {
	    nFwdGen[isel][ieta] = 0;
	    nBwdGen[isel][ieta] = 0;
	  }						      
	}                                                     

	//-- loop over HF selection
	for(int isel = 0; isel < nHFsel; isel++) { 
	  
	  if(!HFgensel[isel]) continue;	     
	  hselection_moca_gen->Fill(3+isel,weight);

	  //-- loop over generated particles
	  for(int igen = 0; igen < ngenparticle; ++igen) {

	    if(!isGoodParticle.at(igen)) continue;
	  
	    //-- loop nfwd - nbwd
	    for(int ieta = 0; ieta < n_delta_eta; ++ieta) {
	      
	      if(ParticleEta->at(igen) > eta_plus_min[ieta] && ParticleEta->at(igen) < eta_plus_max[ieta]) 
		nFwdGen[isel][ieta]+=1;          
	      
	      if(ParticleEta->at(igen) > eta_minus_min[ieta] && ParticleEta->at(igen) < eta_minus_max[ieta]) 
		nBwdGen[isel][ieta]+=1;
	             
	    } //-- end loop nfwd - nbwd
	  } //-- end loop over generated particles
	} //-- end loop over HF selection 

	//-- correlation coefficient
	for(int isel = 0; isel < nHFsel; isel++) {	 
	  for(int ieta = 0; ieta < n_delta_eta; ++ieta) {
	    hcorrelationGen[isel][ieta]->Fill(nFwdGen[isel][ieta],nBwdGen[isel][ieta],weight);	 
	    hcorrelationTheo[isel][ieta]->Fill(nFwdGen[isel][ieta],nBwdGen[isel][ieta],1);
	  } 
	}

      } //-- end loop moca gen 

      //-- HF tower selection
      
      //-- cut 0: down - 1: nominal - 2: up
      bool activity_reco_plus[nHFcut];
      bool activity_reco_minus[nHFcut];

      for(int icut = 0; icut < nHFcut; ++icut) {		
	activity_reco_plus[icut] = false; 				
	activity_reco_minus[icut] = false;	
      }

      //-- cut 0: down - 1: nominal - 2: up
      //-- selection 0: inclusive - 1: inelastic - 2: NSD - 3: SD	
      bool HFrecosel[nHFcut][nHFsel]; 

      for(int icut = 0; icut < nHFcut; ++icut) {		
	for(int isel = 0; isel < nHFsel; isel++) {
	  HFrecosel[icut][isel] = false;   
	}
      }                                               

      //-- different cuts
      for(int icut = 0; icut < nHFcut; ++icut) {
	  
	for(int itower = 0; itower < nHFTower; ++itower) {
	    
	  double EHFtower = HFscaling*HFTowerEnergy->at(itower);
	  double EtaHFtower = HFTowerEta->at(itower);

	  if(CheckHFsel && itower == nHFTower-1) cout<<endl;

	  //-- activity
	  if(!(GetHFrange(EtaHFtower) && EHFtower > EHFcut[icut])) continue;

	  //-- HF plus
	  if(EtaHFtower > 0) {
	    activity_reco_plus[icut] = true;		  
	    if(CheckHFsel) cout<<"HF selection "<<HFcut[icut]<<" - HF cut = "<<EHFcut[icut]<<" GeV"<<endl;
	  }

	  //-- HF minus
	  if(EtaHFtower < 0) {		
	    activity_reco_minus[icut] = true;	    
	    if(CheckHFsel) cout<<"HF selection "<<HFcut[icut]<<" - HF cut = "<<EHFcut[icut]<<" GeV"<<endl;
	  }
 
	} //-- end loop over HFtower
      } //-- end loop over nHFcut

      for(int icut = 0; icut < nHFcut; ++icut) {
	//-- isel = 0 - inclusive: no HF
	HFrecosel[icut][0] = true;   	      
      
	//-- isel = 1 - inelastic: HF OR						   
	if(activity_reco_plus[icut] || activity_reco_minus[icut]) HFrecosel[icut][1] = true;	   
      
	//-- isel = 2 - NSD: HF AND		
	if(activity_reco_plus[icut] && activity_reco_minus[icut]) HFrecosel[icut][2] = true;    

					   
	//-- isel = 3 - SD: HF XOR						
	if(activity_reco_plus[icut] && !activity_reco_minus[icut]) HFrecosel[icut][3] = true;	   
	if(!activity_reco_plus[icut] && activity_reco_minus[icut]) HFrecosel[icut][3] = true;
      }

      if(CheckHFsel) {		
	
	for(int icut = 0; icut < nHFcut; ++icut) {  	 
	  cout<<HFcut[icut]<<endl<<endl;
	  cout<<"activity in HF plus = "<<activity_reco_plus[icut]<<endl;	    
	  cout<<"activity in HF minus = "<<activity_reco_minus[icut]<<endl<<endl;  

	  for(int isel = 0; isel < nHFsel; isel++) {
	    cout<<HFsel[isel]<<" = "<<HFrecosel[icut][isel]<<endl;	  
	  }
	  cout<<endl;
	}
	
	getchar();							   
      }                                                                  
      
      //-- Apply filter data - moca reco - at least one good vertex - vertices, pileup and systematics
      if((isData && data_sel_at_least_one_gv) || (!isData && moca_reco_sel_at_least_one_gv)) {
	
	hVNum->Fill(nVertex,weight);
	hVNumGood->Fill(nGoodVertex15,weight);

	double nFwdPU[nHFsel][n_delta_eta];		      
	double nBwdPU[nHFsel][n_delta_eta];		      
								      
	for(int isel = 0; isel < nHFsel; isel++) {	      
	  for(int ieta = 0; ieta < n_delta_eta; ++ieta) {     
	    nFwdPU[isel][ieta] = 0;			      
	    nBwdPU[isel][ieta] = 0;			      
	  }						      
	}                                                     

	//-- loop over HF selection
	for(int isel = 0; isel < nHFsel; isel++) {  
		
	  //-- cut 1 nominal
	  if(!HFrecosel[1][isel]) continue;            

	  //-- loop on tracks associated to the good vertex with the highest scalar pT sum
	  for (int itrack = 0; itrack < VnTrack->at(iGoodVertex15_max); itrack++) {	      

	    //-- cut on track pT
	    if(VTrackPt->at(iGoodVertex15_max).at(itrack) < 0.5) continue; 

	    //-- good track selection
	    if(!isGoodTrack.at(iGoodVertex15_max).at(itrack)) continue;

	    //-- loop nfwd - nbwd
	    for(int ieta = 0; ieta < n_delta_eta; ++ieta) {
	      
	      if(VTrackEta->at(iGoodVertex15_max).at(itrack) > eta_plus_min[ieta] && VTrackEta->at(iGoodVertex15_max).at(itrack) < eta_plus_max[ieta])	 
		nFwdPU[isel][ieta]+=1;	    
	      
	      if(VTrackEta->at(iGoodVertex15_max).at(itrack) > eta_minus_min[ieta] && VTrackEta->at(iGoodVertex15_max).at(itrack) < eta_minus_max[ieta]) 
		nBwdPU[isel][ieta]+=1;
	      
	    } //-- end loop nfwd - nbwd	     
	  } //-- end loop on tracks associated to the good vertex with the highest scalar pT sum
	} //-- end loop over HF selection	  	       
	
	//-- correlation coefficient						       
	for(int isel = 0; isel < nHFsel; isel++) {	   			       
	  for(int ieta = 0; ieta < n_delta_eta; ++ieta) {			       
	    hcorrelationPU[isel][ieta]->Fill(nFwdPU[isel][ieta],nBwdPU[isel][ieta],weight);  
	  }									       
	}                                                                              

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

	  hPUintVNum->Fill(1+PUint,nGoodVertex15,weight);
	  hPUearlyVNum->Fill(1+PUearly,nGoodVertex15,weight);
	  hPUlateVNum->Fill(1+PUlate,nGoodVertex15,weight);
	  hPUtotVNum->Fill(1+PUtot,nGoodVertex15,weight);     
	}                                           

	if(CheckScalarPtSum) {

	  for (int ivertex = 0; ivertex < nVertex; ivertex++) {
	    
	    //-- filter on good vertex
	    if(!isGoodVertex15[ivertex]) continue;	  
	    cout<<"good vertex 15 (ivertex = "<<ivertex<<"): scalar pT sum = "<<ScalarPtSum[ivertex]<<" GeV"<<endl;
	    
	    //-- filter on good vertex with highest scalar pT sum
	    if(ivertex != iGoodVertex15_max) continue;
	    cout<<"good vertex 15 (ivertex = "<<ivertex<<"): highest scalar pT sum = "<<ScalarPtSum[ivertex]<<" GeV"<<endl;
	  } 
	  
	  getchar();
	  cout<<endl;
	} //-- end CheckScalarPtSum

      } //-- end loop on data - moca reco at least one good vertex

      //-- Apply filter data - moca reco - exactly one good vertex
      if((isData && data_sel) || (!isData && moca_reco_sel)) {

	//-- HF selection
	for(int isel = 0; isel < nHFsel; isel++) {			      
	  //-- cut 1 nominal
	  if(isData && HFrecosel[1][isel]) hselection_data->Fill(7+isel,weight);   
	  //-- cut 1 nominal
	  if(!isData && HFrecosel[1][isel]) hselection_moca_reco->Fill(4+isel,weight);
	}                                                                       

	//-- HF tower 
	int nHFTower_HFcut[nHFcut];
	for(int icut = 0; icut < nHFcut; ++icut) nHFTower_HFcut[icut] = 0;
	
	double HFTowerEtot_HFcut[nHFcut];
	for(int icut = 0; icut < nHFcut; ++icut) HFTowerEtot_HFcut[icut] = 0;
	
	double HFTowerLeadingEnergy_HFcut[nHFcut];
	for(int icut = 0; icut < nHFcut; ++icut) HFTowerLeadingEnergy_HFcut[icut] = 0;

	//-- HF plus tower 								  
	int nHFplusTower_HFcut[nHFcut];						  
	for(int icut = 0; icut < nHFcut; ++icut) nHFplusTower_HFcut[icut] = 0;		  
	
	double HFplusTowerEtot_HFcut[nHFcut];					  
	for(int icut = 0; icut < nHFcut; ++icut) HFplusTowerEtot_HFcut[icut] = 0;	  
	
	double HFplusTowerLeadingEnergy_HFcut[nHFcut];				  
	for(int icut = 0; icut < nHFcut; ++icut) HFplusTowerLeadingEnergy_HFcut[icut] = 0;	

	//-- HF minus tower 								
	int nHFminusTower_HFcut[nHFcut];						  	
	for(int icut = 0; icut < nHFcut; ++icut) nHFminusTower_HFcut[icut] = 0;		
	
	double HFminusTowerEtot_HFcut[nHFcut];					  	
	for(int icut = 0; icut < nHFcut; ++icut) HFminusTowerEtot_HFcut[icut] = 0;	  	
	
	double HFminusTowerLeadingEnergy_HFcut[nHFcut];				  	
	for(int icut = 0; icut < nHFcut; ++icut) HFminusTowerLeadingEnergy_HFcut[icut] = 0;	

	//-- different cuts
	for(int icut = 0; icut < nHFcut; ++icut) {
	  
	  for(int itower = 0; itower < nHFTower; ++itower) {
	    
	    double EHFtower = HFscaling*HFTowerEnergy->at(itower);
	    
	    if(EHFtower > EHFcut[icut]) {
	      
	      //-- HF plus and HF minus together
	      hHFTowerEnergy[icut]->Fill(EHFtower,weight);
	      hHFTowerEta[icut]->Fill(HFTowerEta->at(itower),weight);
	      hHFTowerPhi[icut]->Fill(HFTowerPhi->at(itower),weight);
	      hHFTowerEtaPhi[icut]->Fill(HFTowerEta->at(itower),HFTowerPhi->at(itower),weight);
	      
	      nHFTower_HFcut[icut]++;
	      HFTowerEtot_HFcut[icut]+=EHFtower;
	      
	      if(EHFtower >  HFTowerLeadingEnergy_HFcut[icut]) 
		HFTowerLeadingEnergy_HFcut[icut] = EHFtower;
	      
	      //-- HF plus alone
	      if(HFTowerEta->at(itower) > 0) {		
		 
		hHFplusTowerEnergy[icut]->Fill(EHFtower,weight);
		hHFplusTowerEta[icut]->Fill(HFTowerEta->at(itower),weight);
		hHFplusTowerPhi[icut]->Fill(HFTowerPhi->at(itower),weight);
		hHFplusTowerEtaPhi[icut]->Fill(HFTowerEta->at(itower),HFTowerPhi->at(itower),weight);

		nHFplusTower_HFcut[icut]++;
		HFplusTowerEtot_HFcut[icut]+=EHFtower;
	      
		if(EHFtower >  HFplusTowerLeadingEnergy_HFcut[icut]) 
		  HFplusTowerLeadingEnergy_HFcut[icut] = EHFtower;     
	      }

	      //-- HF minus alone
	      if(HFTowerEta->at(itower) < 0) {		

		hHFminusTowerEnergy[icut]->Fill(EHFtower,weight);
		hHFminusTowerEta[icut]->Fill(HFTowerEta->at(itower),weight);
		hHFminusTowerPhi[icut]->Fill(HFTowerPhi->at(itower),weight);
		hHFminusTowerEtaPhi[icut]->Fill(HFTowerEta->at(itower),HFTowerPhi->at(itower),weight);

		nHFminusTower_HFcut[icut]++;
		HFminusTowerEtot_HFcut[icut]+=EHFtower;
		
		if(EHFtower >  HFminusTowerLeadingEnergy_HFcut[icut]) 
		  HFminusTowerLeadingEnergy_HFcut[icut] = EHFtower;
	      }
 
	    } //-- end EHFcut
	  } //-- end loop over HFtower

	  //-- HF plus and HF minus together
	  hnHFTower[icut]->Fill(nHFTower_HFcut[icut],weight);
	  hHFTowerEtot[icut]->Fill(HFTowerEtot_HFcut[icut],weight);
	  if(HFTowerLeadingEnergy_HFcut[icut] > 0) hHFTowerLeadingEnergy[icut]->Fill(HFTowerLeadingEnergy_HFcut[icut],weight);

	  //-- HF plus alone
	  hnHFplusTower[icut]->Fill(nHFplusTower_HFcut[icut],weight);			   
	  hHFplusTowerEtot[icut]->Fill(HFplusTowerEtot_HFcut[icut],weight);		   
	  if(HFplusTowerLeadingEnergy_HFcut[icut] > 0) hHFplusTowerLeadingEnergy[icut]->Fill(HFplusTowerLeadingEnergy_HFcut[icut],weight);

	  //-- HF minus alone								   
	  hnHFminusTower[icut]->Fill(nHFminusTower_HFcut[icut],weight);			   
	  hHFminusTowerEtot[icut]->Fill(HFminusTowerEtot_HFcut[icut],weight);		   
	  if(HFminusTowerLeadingEnergy_HFcut[icut] > 0) hHFminusTowerLeadingEnergy[icut]->Fill(HFminusTowerLeadingEnergy_HFcut[icut],weight);
	
	} //-- end loop over nHFcut
	
	//-- vertex and tracks
	
	double nFwd[nHFsel][n_delta_eta];
	double nBwd[nHFsel][n_delta_eta];

	double nFwdDown[nHFsel][n_delta_eta];
	double nBwdDown[nHFsel][n_delta_eta];

	double nFwdUp[nHFsel][n_delta_eta];
	double nBwdUp[nHFsel][n_delta_eta];

	double nFwdEffi[nHFsel][n_delta_eta];		      
	double nBwdEffi[nHFsel][n_delta_eta];		      
	
	for(int isel = 0; isel < nHFsel; isel++) {
	  for(int ieta = 0; ieta < n_delta_eta; ++ieta) {
	    nFwd[isel][ieta] = 0;
	    nBwd[isel][ieta] = 0;

	    nFwdDown[isel][ieta] = 0;
	    nBwdDown[isel][ieta] = 0;

	    nFwdUp[isel][ieta] = 0;
	    nBwdUp[isel][ieta] = 0;

	    nFwdEffi[isel][ieta] = 0;			      
	    nBwdEffi[isel][ieta] = 0;			      
	  }						      
	}                                                     
       
	double leading_track_pT = 0;
	double leading_track_eta = 0;
	double leading_track_phi = 0;

	//-- exactly one good vertex
	hBSx->Fill(BSx,weight);
	hBSy->Fill(BSy,weight);
	hBSz->Fill(BSz,weight);
	
	nSelTrack100 = 0;
	nSelTrack96 = 0;  

	//-- loop on tracks - exactly one good vertex
	for (int itrack = 0; itrack < VnTrack->at(0); itrack++) {	      

	  double trackeffi = random->Rndm();

	  //-- cut on track pT
	  if(VTrackPt->at(0).at(itrack) < 0.5) continue; 

	  //-- fill quality before selection 
	  hTrackQuality->Fill(VTrackQuality->at(0).at(itrack),weight);

	  //-- good track selection
	  if(!isGoodTrack.at(0).at(itrack)) continue;

	  nSelTrack100+=1;
	  if(trackeffi < trackeffi_cut) nSelTrack96+=1;

	  //-- loop over HF selection
	  for(int isel = 0; isel < nHFsel; isel++) {

	    //-- cut 0 down													
	    if(HFrecosel[0][isel]) {												
																
	      //-- loop over nfw - nbwd												
	      for(int ieta = 0; ieta < n_delta_eta; ++ieta) {									
																
		if(VTrackEta->at(0).at(itrack) > eta_plus_min[ieta] && VTrackEta->at(0).at(itrack) < eta_plus_max[ieta])	
		  nFwdDown[isel][ieta]+=1;          										
  																
		if(VTrackEta->at(0).at(itrack) > eta_minus_min[ieta] && VTrackEta->at(0).at(itrack) < eta_minus_max[ieta]) 	
		  nBwdDown[isel][ieta]+=1;												
	      															
	      } //-- end loop over nfwd - nbwd											
	    } //-- end cut 0 down                                                                                            


	    //-- cut 1 nominal
	    if(HFrecosel[1][isel]) {

	      //-- loop over nfw - nbwd
	      for(int ieta = 0; ieta < n_delta_eta; ++ieta) {

		if(VTrackEta->at(0).at(itrack) > eta_plus_min[ieta] && VTrackEta->at(0).at(itrack) < eta_plus_max[ieta]) {		
		  nFwd[isel][ieta]+=1;          
		  if(trackeffi < trackeffi_cut) nFwdEffi[isel][ieta]+=1;
		}
  
		if(VTrackEta->at(0).at(itrack) > eta_minus_min[ieta] && VTrackEta->at(0).at(itrack) < eta_minus_max[ieta]) {
		  nBwd[isel][ieta]+=1;
		  if(trackeffi < trackeffi_cut) nBwdEffi[isel][ieta]+=1;
		}
	      
	      } //-- end loop over nfwd - nbwd
	    } //-- end cut 1 nominal                                                                                            


	    //-- cut 2 up
	    if(HFrecosel[2][isel]) {												
																
	      //-- loop over nfw - nbwd												
	      for(int ieta = 0; ieta < n_delta_eta; ++ieta) {									
																
		if(VTrackEta->at(0).at(itrack) > eta_plus_min[ieta] && VTrackEta->at(0).at(itrack) < eta_plus_max[ieta]) 
		  nFwdUp[isel][ieta]+=1;          										
  																
		if(VTrackEta->at(0).at(itrack) > eta_minus_min[ieta] && VTrackEta->at(0).at(itrack) < eta_minus_max[ieta]) 	
		  nBwdUp[isel][ieta]+=1;												
	      															
	      } //-- end loop over nfwd - nbwd											
	    } //-- end cut 2 up


	  } //-- end loop over HF selection
	 	       
	  hTrackPt->Fill(VTrackPt->at(0).at(itrack),weight);
	  hTrackEta->Fill(VTrackEta->at(0).at(itrack),weight);
	  hTrackPhi->Fill(VTrackPhi->at(0).at(itrack),weight);

	  if(TMath::Abs(VTrackEta->at(0).at(itrack)) < 1) hTrackPhiCentral->Fill(VTrackPhi->at(0).at(itrack),weight);
	  if(TMath::Abs(VTrackEta->at(0).at(itrack)) > 1) hTrackPhiFwd->Fill(VTrackPhi->at(0).at(itrack),weight);

	  hTrackEtaPhi->Fill(VTrackEta->at(0).at(itrack),VTrackPhi->at(0).at(itrack),weight);

	  if(VTrackPt->at(0).at(itrack) > leading_track_pT) {
	    leading_track_pT = VTrackPt->at(0).at(itrack);
	    leading_track_eta = VTrackEta->at(0).at(itrack);
	    leading_track_phi = VTrackPhi->at(0).at(itrack);
	  }
	  
	  double abseta = TMath::Abs(VTrackEta->at(0).at(itrack));

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
	    hTrackdxyNum2D->Fill(nGoodTrack[0],dxysigmaxy,weight); 
	  }
	  
	  if(VTrackdzError->at(0).at(itrack) != 0) {
	    double dzsigmaz = VTrackdz->at(0).at(itrack)/VTrackdzError->at(0).at(itrack);
	    hTrackdz->Fill(TMath::Abs(dzsigmaz),weight);
	    hTrackdzEta->Fill(VTrackEta->at(0).at(itrack),TMath::Abs(dzsigmaz),weight);
	    hTrackdzNhit->Fill(VTrackValidHit->at(0).at(itrack),TMath::Abs(dzsigmaz),weight);
	    
	    hTrackdzEta2D->Fill(VTrackEta->at(0).at(itrack),dzsigmaz,weight); 
	    hTrackdzNum2D->Fill(nGoodTrack[0],dzsigmaz,weight); 
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
	  
	} //-- end loop on tracks - exactly one good vertex

	hTrackNum->Fill(nGoodTrack[0],weight);

	//-- exactly one good vertex
	hVx->Fill(VxBS[0],weight);
	hVy->Fill(VyBS[0],weight);
	hVxy->Fill(VxBS[0],VyBS[0],weight);
	hVz->Fill(Vz->at(0),weight);
	hVrho->Fill(rhoBS[0],weight);

	hVchi2->Fill(Vchi2->at(0),weight);
	hVndof->Fill(Vndof->at(0),weight);
	hVchi2ndof->Fill(Vchi2ndof->at(0),weight);

	hVzError->Fill(VzError->at(0),weight);
	hVzErrorVz->Fill(Vz->at(0),VzError->at(0),weight);
	hVzErrorTrackNum->Fill(nGoodTrack[0],VzError->at(0),weight);
	
	//-- leading track
	if(leading_track_pT > 0) {
	  hLeadingTrackPt->Fill(leading_track_pT,weight);
	  hLeadingTrackEta->Fill(leading_track_eta,weight);
	  hLeadingTrackPhi->Fill(leading_track_phi,weight);
	}
	
	//-- correlation coefficient
	for(int isel = 0; isel < nHFsel; isel++) {	   
	  for(int ieta = 0; ieta < n_delta_eta; ++ieta) {
	    hcorrelation[isel][ieta]->Fill(nFwd[isel][ieta],nBwd[isel][ieta],weight);
	    hcorrelationDown[isel][ieta]->Fill(nFwdDown[isel][ieta],nBwdDown[isel][ieta],weight);
	    hcorrelationUp[isel][ieta]->Fill(nFwdUp[isel][ieta],nBwdUp[isel][ieta],weight);
	    hcorrelationEffi[isel][ieta]->Fill(nFwdEffi[isel][ieta],nBwdEffi[isel][ieta],weight);
	  }
	}                                                                              

	if(CheckEffi) {

	  for(int i = 0; i < n_delta_eta; ++i) {
	    cout<<"bin "<<i+1<<endl;
	    cout<<"nFwd = "<<nFwd[i]<<" - nFwd effi sys = "<<nFwdEffi[i]<<endl;
	    cout<<"nBwd = "<<nBwd[i]<<" - nBwd effi sys = "<<nBwdEffi[i]<<endl;
	  }
	  
	  cout<<endl;
	  cout<<"n selected tracks 100% efficiency = "<<nSelTrack100<<endl;
	  cout<<"n selected tracks 96% efficiency = "<<nSelTrack96<<endl;
	  if(nSelTrack100 > 0) cout<<"ratio = "<<100*nSelTrack96/nSelTrack100<<" %"<<endl<<endl;
	}

      } //-- end loop on data - moca reco exactly one good vertex

      isGoodTrack.clear();
      if(!isData) isGoodParticle.clear();

      nSelTrack100_tot+=nSelTrack100;
      nSelTrack96_tot+=nSelTrack96;

    } //-- end event loop
    
    delete tree;
    file->Close();
    
  } //-- end file loop
  
  //-- HF tower eta - phi normalization 
  for(int icut = 0; icut < nHFcut; ++icut) {
    double scaling = 0;
    if(!isData) scaling = 1/Nmoca_reco_sel;
    if(isData) scaling = 1/Ndata_sel;

    hHFTowerEtaPhi[icut]->Scale(scaling);
    hHFplusTowerEtaPhi[icut]->Scale(scaling);
    hHFminusTowerEtaPhi[icut]->Scale(scaling);
  }

  //-- correlation coefficient
  GetCorrelation(hcorrelation,hrho);

  GetCorrelation(hcorrelationDown,hrhoDown);
  GetCorrelation(hcorrelationUp,hrhoUp);

  GetCorrelation(hcorrelationPU,hrhoPU);
  GetCorrelation(hcorrelationEffi,hrhoEffi);
			 
  if(!isData) GetCorrelation(hcorrelationGen,hrhoGen);
  if(!isData) GetCorrelation(hcorrelationTheo,hrhoTheo);

  if(!isData) ComputeHistoRatio(hrhoGen,hrho,hrhoCF);
  if(!isData) ComputeHistoRatio(hrhoGen,hrhoPU,hrhoCFPU);
  
  cout<<endl;
  
  //-- Check histo 
  CheckHistoSelection();
  
  CheckHistoHFTower();
  CheckHistoHFplusTower();
  CheckHistoHFminusTower();
  
  CheckHistoTrack();
  CheckHistoVertex();
  CheckHistoBS();

  CheckHistoPU(isData);
  CheckHistoRho(isData);

  //-- write histo to output root file 
  Char_t outputfile_name[200];
  sprintf(outputfile_name,"%s",file_name.Data());
  TFile* foutput = new TFile(outputfile_name,"RECREATE");
  foutput->cd();

  //-- Write histo
  WriteHistoSelection();
  
  WriteHistoHFTower();
  WriteHistoHFplusTower();
  WriteHistoHFminusTower();

  WriteHistoTrack();
  WriteHistoVertex();
  WriteHistoBS();

  WriteHistoPU(isData);
  WriteHistoRho(isData);
					        
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
    cout<<"number of events inclusive selection: "<<hselection_data->GetBinContent(7)<<endl;
    cout<<"number of events inelastic selection: "<<hselection_data->GetBinContent(8)<<endl;
    cout<<"number of events NSD selection: "<<hselection_data->GetBinContent(9)<<endl;
    cout<<"number of events SD selection: "<<hselection_data->GetBinContent(10)<<endl;       
  }
   
  if(!isData) {
    cout<<endl<<"moca at reco level: "<<endl;
    cout<<"number of events before selection: "<<hselection_moca_reco->GetBinContent(1)<<endl;
    cout<<"number of events after at least one good vertex: "<<hselection_moca_reco->GetBinContent(2)<<endl;
    cout<<"number of events after exactly one good vertex: "<<hselection_moca_reco->GetBinContent(3)<<endl;
    cout<<"number of events inclusive selection: "<<hselection_moca_reco->GetBinContent(4)<<endl; 
    cout<<"number of events inelastic selection: "<<hselection_moca_reco->GetBinContent(5)<<endl; 
    cout<<"number of events NSD selection: "<<hselection_moca_reco->GetBinContent(6)<<endl;	     
    cout<<"number of events SD selection: "<<hselection_moca_reco->GetBinContent(7)<<endl;       

    cout<<endl<<"moca at gen level: "<<endl;
    cout<<"number of events before selection: "<<hselection_moca_gen->GetBinContent(1)<<endl;
    cout<<"number of events after vertex selection: "<<hselection_moca_gen->GetBinContent(2)<<endl;
    cout<<"number of events inclusive selection: "<<hselection_moca_gen->GetBinContent(3)<<endl;
    cout<<"number of events inelastic selection: "<<hselection_moca_gen->GetBinContent(4)<<endl;
    cout<<"number of events NSD selection: "<<hselection_moca_gen->GetBinContent(5)<<endl;
    cout<<"number of events SD selection: "<<hselection_moca_gen->GetBinContent(6)<<endl;               
  }
  
  cout<<endl;
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
  if(!isData) cout<<"selected events in mc at reco level - only pu weight: "<<Nmoca_reco_sel_no_vz_weight<<endl;
  if(!isData) cout<<"selected events in mc at gen level - global weight: "<<Nmoca_gen_sel<<endl<<endl;

  long double ratio_moca_reco_sel_1 = 0;						       
  if(!isData) ratio_moca_reco_sel_1 = Nmoca_reco_sel_no_weight/Nmoca_reco_sel;		   	       
  if(!isData) cout<<"no weight / global weight: "<<ratio_moca_reco_sel_1<<endl;              
  
  long double ratio_moca_reco_sel_2 = 0;						       
  if(!isData) ratio_moca_reco_sel_2 = Nmoca_reco_sel_no_weight/Nmoca_reco_sel_no_pu_weight;  	       
  if(!isData) cout<<"no weight / only vertex weight: "<<ratio_moca_reco_sel_2<<endl;         
  
  long double ratio_moca_reco_sel_3 = 0;						       
  if(!isData) ratio_moca_reco_sel_3 = Nmoca_reco_sel_no_weight/Nmoca_reco_sel_no_vz_weight;  	       
  if(!isData) cout<<"no weight / only pu weight: "<<ratio_moca_reco_sel_3<<endl<<endl;       
  
  cout<<"n selected tracks 100% efficiency = "<<nSelTrack100_tot<<endl;
  cout<<"n selected tracks 96% efficiency = "<<nSelTrack96_tot<<endl;
  if(nSelTrack100_tot > 0) cout<<"ratio = "<<100*nSelTrack96_tot/nSelTrack100_tot<<" %"<<endl<<endl;

  cout<<"correlation coefficient nominal"<<endl;

  for(int isel = 0; isel < nHFsel; isel++) {

    cout<<endl;
    cout<<HFsel[isel]<<endl<<endl;
    
    for (int i = 1; i <= hrho[isel]->GetNbinsX(); i++) {
      cout<<"rho bin "<<i<<": "<<hrho[isel]->GetBinContent(i)<<" +/- "<<hrho[isel]->GetBinError(i)<<endl;
    }
  
  }                  
 
  if(!isData) {													      					       
    cout<<endl;
    cout<<"correction factor"<<endl;

    for(int isel = 0; isel < nHFsel; isel++) {									       

      cout<<endl;
      cout<<"correction factor - "<<HFsel[isel]<<endl<<endl;                             								   

      for (int i = 1; i <= hrhoCF[isel]->GetNbinsX(); i++) { 						       
  	cout<<"CF bin "<<i<<": "<<hrhoCF[isel]->GetBinContent(i)<<" +/- "<<hrhoCF[isel]->GetBinError(i)<<endl;	             
      }							
    }							       
    														       
    cout<<endl;													       
    cout<<"correction factor PU"<<endl;

    for(int isel = 0; isel < nHFsel; isel++) {									       

      cout<<endl;
      cout<<"correction factor PU - "<<HFsel[isel]<<endl<<endl;                             								 

      for (int i = 1; i <= hrhoCFPU[isel]->GetNbinsX(); i++) {						       
        cout<<"CF PU bin "<<i<<": "<<hrhoCFPU[isel]->GetBinContent(i)<<" +/- "<<hrhoCFPU[isel]->GetBinError(i)<<endl;    
      }							
    }							       

    cout<<endl;
  }
    
  TString temp = outputfile_name;
  TString logfile_name(temp.Remove(temp.Length()-4));
  logfile_name+="txt";
  
  //-- logfile
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
    logfile<<"number of events inclusive selection: "<<hselection_data->GetBinContent(7)<<endl; 
    logfile<<"number of events inelastic selection: "<<hselection_data->GetBinContent(8)<<endl; 
    logfile<<"number of events NSD selection: "<<hselection_data->GetBinContent(9)<<endl;	     
    logfile<<"number of events SD selection: "<<hselection_data->GetBinContent(10)<<endl;       
  }
  
  if(!isData) {
    logfile<<endl<<"moca at reco level: "<<endl;
    logfile<<"number of events before selection: "<<hselection_moca_reco->GetBinContent(1)<<endl;
    logfile<<"number of events after at least one good vertex: "<<hselection_moca_reco->GetBinContent(2)<<endl;
    logfile<<"number of events after exactly one good vertex: "<<hselection_moca_reco->GetBinContent(3)<<endl;
    logfile<<"number of events inclusive selection: "<<hselection_moca_reco->GetBinContent(4)<<endl;
    logfile<<"number of events inelastic selection: "<<hselection_moca_reco->GetBinContent(5)<<endl;
    logfile<<"number of events NSD selection: "<<hselection_moca_reco->GetBinContent(6)<<endl;	 
    logfile<<"number of events SD selection: "<<hselection_moca_reco->GetBinContent(7)<<endl;       

    logfile<<endl<<"moca at gen level: "<<endl;								
    logfile<<"number of events before selection: "<<hselection_moca_gen->GetBinContent(1)<<endl;		
    logfile<<"number of events after vertex selection: "<<hselection_moca_gen->GetBinContent(2)<<endl;	
    logfile<<"number of events inclusive selection: "<<hselection_moca_gen->GetBinContent(3)<<endl;	
    logfile<<"number of events inelastic selection: "<<hselection_moca_gen->GetBinContent(4)<<endl;	
    logfile<<"number of events NSD selection: "<<hselection_moca_gen->GetBinContent(5)<<endl;		
    logfile<<"number of events SD selection: "<<hselection_moca_gen->GetBinContent(6)<<endl;               
  }
  
  logfile<<endl;
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
  if(!isData) logfile<<"selected events in mc at reco level - only pu weight: "<<Nmoca_reco_sel_no_vz_weight<<endl;
  if(!isData) logfile<<"selected events in mc at gen level - global weight: "<<Nmoca_gen_sel<<endl<<endl;  

  if(!isData) logfile<<"no weight / global weight: "<<ratio_moca_reco_sel_1<<endl;
  if(!isData) logfile<<"no weight / only vertex weight: "<<ratio_moca_reco_sel_2<<endl;
  if(!isData) logfile<<"no weight / only pu weight: "<<ratio_moca_reco_sel_3<<endl<<endl;
  
  logfile.close();
  cout<<endl<<"log file: "<<logfile_name<<" created"<<endl;
}



