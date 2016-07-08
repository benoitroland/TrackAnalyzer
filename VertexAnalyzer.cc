#include "VertexAnalyzer.h"

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

#define BScheck 0
#define VertexDebug 0

VertexAnalyzer::VertexAnalyzer() { }

VertexAnalyzer::~VertexAnalyzer() { }

void VertexAnalyzer::Loop(TString inputdir,TObjArray* filelist,TString type, int iteration, TString file_name) {

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

  cout<<endl;
  cout<<"type = "<<type<<endl;
  cout<<"iteration = "<<iteration<<endl<<endl;
  cout<<"press enter to continue"<<endl;
  getchar();
  
  //-- event
  long double Nevt_tot = 0;
  long double Nevt_tot_no_weight = 0;

  long double Ndata_sel = 0;
  long double Ndata_sel_no_weight = 0;

  long double Nmoca_reco_sel = 0;
  long double Nmoca_reco_sel_no_weight = 0;

  //-- histo label and title

  TString label,title;
  TString xleg, yleg;

  //-- selection
  TH1D *hselection_data = MakeHisto("hselection_data","Event Selection Data","step","N evts",6,0.5,6.5);
  TH1D *hselection_moca_reco = MakeHisto("hselection_moca_reco","Event Selection MC reco","step","N evts",3,0.5,3.5);
  cout<<endl;

  //-- vertex distributions exactly one good vertex
  TH1D *hVx = MakeHisto("hVx","vertex x","vertex x [cm]","N evts",200,-0.21,0.19);
  TH1D *hVy = MakeHisto("hVy","vertex y","vertex y [cm]","N evts",200,-0.21,0.19);
  TH2D *hVxy = MakeHisto("hVxy","vertex y versus x","vertex x [cm]","vertex y [cm]",200,-0.21,0.19,200,-0.21,0.19);
  TH1D *hVz = MakeHisto("hVz","vertex z","vertex z [cm]","N evts",210,-21.5,20.5);
  TH1D *hVrho = MakeHisto("hVrho","vertex rho","vertex #rho [cm]","N evts",2000,0,1);
  TH1D *hVrhoBS = MakeHisto("hVrhoBS","vertex rho BS","vertex #rho BS [cm]","N evts",10000,0,1);
  TH1D *hVchi2 = MakeHisto("hVchi2","vertex chi2","vertex #chi^{2}","N evts",500,0,100);
  TH1D *hVndof = MakeHisto("hVndof","vertex ndof","vertex ndof","N evts",100,0,100);
  TH1D *hVchi2ndof = MakeHisto("hVchi2ndof","vertex chi2/ndof","vertex #chi^{2}/ndof","N evts",200,0,10); 
  TH1D *hVzError = MakeHisto("hVzError","vertex z error","vertex z error [cm]","N evts",200,0,0.2);
  TH2D *hVzErrorVz = MakeHisto("hVzErrorVz","vertex z error versus vertex z","vertex z [cm]","vertex z error [cm]",100,-15,15,100,0,0.2);
  cout<<endl;

  //-- beam spot distributions
  TH1D *hBSx = MakeHisto("hBSx","BS x","BS x [cm]","N evts",1000,-0.21,0.19);
  TH1D *hBSy = MakeHisto("hBSy","BS y","BS y [cm]","N evts",1000,-0.21,0.19);
  TH1D *hBSz = MakeHisto("hBSz","BS z","BS z [cm]","N evts",1000,-20.05,19.95);

  //-- weight  	    
  double weight = 1;

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
      cout<<"Error in VertexAnalyzer: could not open file "<<temp_itfile->GetString()<<endl;
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
      cout<<"Error in VertexAnalyzer: could not open file "<<itfile->GetString()<<endl; 
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

      for (int ivert = 0; ivert < nVertex; ++ivert) {	
	double deltax = Vx->at(ivert) - BSx;
	double deltay = Vy->at(ivert) - BSy;
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

      for (int ivertex = 0; ivertex < nVertex; ++ivertex) {	
	
	if(VisFake->at(ivertex)) continue;
	if(!VisValid->at(ivertex)) continue;
	
	if(rhoBS[ivertex] > 0.2) continue;  //-- rhoBS <= 0.2
	if(TMath::Abs(Vz->at(ivertex)) > 20) continue;   //-- |Vz| <= 20

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

	  if(VTrackValidPixelHit->at(ivertex).at(itrack) < 3) continue;  //-- n pixel hits >= 3
	  ngood_track++;
	} //-- end loop over associated tracks
	
	if(ngood_track < 2) continue; //-- nTrack >= 2

	at_least_one_gv = true;
      } //-- end loop over vertex
      
      //-- vertex information - exactly one good vertex
      bool exactly_one_gv = false;					      
      if(at_least_one_gv == true && nVertex == 1) exactly_one_gv = true;

      //-- vertex weight montecarlo 

      weight = 1;

      //-- iteration 0: to determine weight
      //-- iteration 1: to check weight and determine global factor
      //-- iteration 2: to check global factor

      if(!isData && exactly_one_gv) weight = GetVertexWeight(Vz->at(0),AllMC,iteration); 

      if(VertexDebug) {                                                     
	cout<<"n vertex = "<<nVertex<<" - vertex weight = "<<weight<<endl;
	getchar();
      }

      //-- number of events 
      
      Nevt_tot+=weight;
      Nevt_tot_no_weight+=1;

      if(isData) hselection_data->Fill(1,weight);
      if(!isData) hselection_moca_reco->Fill(1,weight);

      if((int (Nevt_tot_no_weight))%10000 == 0) cout<<"total number of events done = "<<Nevt_tot<<" ("<<100*Nevt_tot/Nevt_all_files<<"%)"<<endl;
      if((int (Nevt_tot_no_weight))%10000 == 0 && isData) cout<<"number of selected events in data = "<<Ndata_sel<<endl<<endl;
      if((int (Nevt_tot_no_weight))%10000 == 0 && !isData) cout<<"number of selected events in mc at reco level = "<<Nmoca_reco_sel<<endl<<endl;

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
	}
      } //-- end loop on moca reco

      //-- Apply filter data - moca reco - exactly one good vertex
      if((isData && data_sel) || (!isData && moca_reco_sel)) {
	
	if(isData && LumiSection < 90) cout<<"warning: "<<"LumiSection = "<<LumiSection<<" < 90"<<endl;

	//-- exactly one good vertex
	hBSx->Fill(BSx,weight);
	hBSy->Fill(BSy,weight);
	hBSz->Fill(BSz,weight);

	//-- exactly one good vertex
	hVx->Fill(Vx->at(0),weight);
	hVy->Fill(Vy->at(0),weight);
	hVxy->Fill(Vx->at(0),Vy->at(0),weight);
	hVz->Fill(Vz->at(0),weight);
	hVrho->Fill(Vrho->at(0),weight);
	hVrhoBS->Fill(rhoBS[0],weight);
	hVchi2->Fill(Vchi2->at(0),weight);
	hVndof->Fill(Vndof->at(0),weight);
	hVchi2ndof->Fill(Vchi2ndof->at(0),weight);
	hVzError->Fill(VzError->at(0),weight);
	hVzErrorVz->Fill(Vz->at(0),VzError->at(0),weight);
	
      } //-- end loop on data - moca reco exactly one good vertex

    } //-- end event loop
  
    delete tree;
    file->Close();
  
  } //-- end file loop

  //-- Check Histo

  //-- selection
  cout<<endl;

  CheckHisto(hselection_data);
  CheckHisto(hselection_moca_reco);
  
  //-- beam spot
  cout<<endl;

  CheckHisto(hBSx);
  CheckHisto(hBSy);
  CheckHisto(hBSz);

  //-- vertex distributions exactly one good vertex
  cout<<endl;

  CheckHisto(hVx);
  CheckHisto(hVy);
  CheckHisto(hVxy);
  CheckHisto(hVz);
  CheckHisto(hVrho);
  CheckHisto(hVrhoBS);
  CheckHisto(hVchi2);
  CheckHisto(hVndof);
  CheckHisto(hVchi2ndof);
  CheckHisto(hVzError);
  CheckHisto(hVzErrorVz);

  //-- write histo to output root file 
  
  Char_t outputfile_name[200];
  sprintf(outputfile_name,"%s",file_name.Data());
  TFile* foutput = new TFile(outputfile_name,"RECREATE");
  foutput->cd();
  
  //-- selection		   
  hselection_data->Write();	   
  hselection_moca_reco->Write();

  //-- beam spot   
  hBSx->Write();
  hBSy->Write();
  hBSz->Write();

  //-- vertex distributions exactly one good vertex
  cout<<endl;
  hVx->Write();
  hVy->Write();
  hVxy->Write();
  hVz->Write();
  hVrho->Write();
  hVrhoBS->Write();
  hVchi2->Write();
  hVndof->Write();
  hVchi2ndof->Write();
  hVzError->Write();
  hVzErrorVz->Write();

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

  cout<<"total number of events - vertex weight: "<<Nevt_tot<<endl;
  cout<<"total number of events - no vertex weight: "<<Nevt_tot_no_weight<<endl;
  long double ratio_evt_tot = Nevt_tot_no_weight/Nevt_tot;
  cout<<"no vertex weight / vertex weight: "<<ratio_evt_tot<<endl<<endl; 

  if(isData) cout<<"selected events in data - vertex weight: "<<Ndata_sel<<endl;
  if(isData) cout<<"selected events in data - no vertex weight: "<<Ndata_sel_no_weight<<endl;
  long double ratio_data_sel = 0;
  if(isData) ratio_data_sel = Ndata_sel_no_weight/Ndata_sel;
  if(isData) cout<<"no vertex weight / vertex weight: "<<ratio_data_sel<<endl<<endl;

  if(!isData) cout<<"selected events in mc at reco level - vertex weight: "<<Nmoca_reco_sel<<endl;
  if(!isData) cout<<"selected events in mc at reco level - no vertex weight: "<<Nmoca_reco_sel_no_weight<<endl;
  long double ratio_moca_reco_sel = 0;
  if(!isData) ratio_moca_reco_sel = Nmoca_reco_sel_no_weight/Nmoca_reco_sel;
  if(!isData) cout<<"no vertex weight / vertex weight: "<<ratio_moca_reco_sel<<endl<<endl;
}



