#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include "MainAnalyzer.h"

void CreateDirectory(struct stat dir, TString dirname) {
  using namespace std;
  if(stat(dirname, &dir) == 0 && S_ISDIR(dir.st_mode)) cout<<"directory "<<dirname<<" is already present!"<<endl;
  else {
    if(mkdir(dirname,0755)==0) cout<<"directory "<<dirname<<" created successfully!"<<endl;
  }
  cout<<"press enter to continue"<<endl;
  getchar();
}

int main(int argc, char *argv[]) {

  struct stat dir;

  using namespace std;

  printf("argc = %d\n", argc);
  for (int i = 0; i<argc; i++){
    static int k = 0;
    printf("argv[%d] = %s\n", k, argv[k]);
    k++;
  }

  if(argc < 3) {
    cout<<endl;
    cout<<"usage: "<<"./Run \"analysis\" \"type\" \"selection or iteration\" "<<endl<<endl;
    cout<<"short cut: data = data-38T - p8  = PYTHIA8-CUETP8M1-38T-PU1p3 - hpp = HERWIGPP-CUETHS1-38T-PU1p3 - epos = EPOS-38T-PU1p3"<<endl<<endl;

    cout<<"analysis: track"<<endl<<endl;
    cout<<"type: data-38T - data-HFRereco-38T - data-PromptReco-38T"<<endl;
    cout<<"      PYTHIA8-CUETP8M1-38T-PU1p3 - PYTHIA8-CUETP8M1-38T-noPU"<<endl;
    cout<<"      HERWIGPP-CUETHS1-38T-PU1p3"<<endl;
    cout<<"      EPOS-38T-PU1p3"<<endl<<endl;
    cout<<"selection: nominal - HFup - HFdown - no_pu_reweight - no_pixel_cut - no_vz_reweight - combined"<<endl<<endl;

    cout<<"analysis: correlation"<<endl<<endl;
    cout<<"type: data-38T - data-HFRereco-38T - data-PromptReco-38T"<<endl;
    cout<<"      PYTHIA8-CUETP8M1-38T-PU1p3 - PYTHIA8-CUETP8M1-38T-noPU"<<endl;
    cout<<"      HERWIGPP-CUETHS1-38T-PU1p3"<<endl;
    cout<<"      EPOS-38T-PU1p3"<<endl<<endl;
    cout<<"selection: nominal - HFup - HFdown - no_pu_reweight - no_pixel_cut - no_vz_reweight"<<endl<<endl;

    cout<<"analysis: vertex"<<endl<<endl;
    cout<<"type: data-38T "<<endl;
    cout<<"      PYTHIA8-CUETP8M1-38T-PU1p3"<<endl;
    cout<<"      HERWIGPP-CUETHS1-38T-PU1p3"<<endl;
    cout<<"      EPOS-38T-PU1p3"<<endl<<endl;
    cout<<"iteration: 0 - 1 - 2"<<endl<<endl;              

    cout<<"analysis: pu"<<endl<<endl;	     
    cout<<"type: data-38T "<<endl;		     
    cout<<"      PYTHIA8-CUETP8M1-38T-PU1p3"<<endl;  
    cout<<"      HERWIGPP-CUETHS1-38T-PU1p3"<<endl;  
    cout<<"      EPOS-38T-PU1p3"<<endl<<endl;	     
    cout<<"iteration: 0 - 1 - 2"<<endl<<endl;        

    return(0);
  }

  TString analysis = TString(argv[1]);
  TString type = TString(argv[2]);

  if(strcmp(type,"data") == 0)  type = "data-38T";
  if(strcmp(type,"p8")==0) type = "PYTHIA8-CUETP8M1-38T-PU1p3";
  if(strcmp(type,"hpp")==0) type = "HERWIGPP-CUETHS1-38T-PU1p3";
  if(strcmp(type,"epos")==0) type = "EPOS-38T-PU1p3";

  TString selection = "";
  if(strcmp(analysis,"track") == 0) selection = TString(argv[3]);
  if(strcmp(analysis,"correlation") == 0) selection = TString(argv[3]);
 
  int iteration = 0;
  if(strcmp(analysis,"vertex") == 0) iteration = atoi(argv[3]);
  if(strcmp(analysis,"pu") == 0) iteration = atoi(argv[3]);

  TString dir_input = "/nfs/dust/cms/user/roland/TrackAnalysis/OutputTree/";    
  TString dir_output =  "/nfs/dust/cms/user/roland/TrackAnalysis/OutputHisto";             

  //-- input directory                                                                                                                                                       
  if(strcmp(type,"data-38T") == 0)            dir_input+="data-38T/";
  if(strcmp(type,"data-HFRereco-38T") == 0)   dir_input+="data-HFRereco-38T/";
  if(strcmp(type,"data-PromptReco-38T") == 0) dir_input+="data-PromptReco-38T/";

  if(strcmp(type,"PYTHIA8-CUETP8M1-38T-noPU")==0)  dir_input+="PYTHIA8-CUETP8M1-38T/";
  if(strcmp(type,"PYTHIA8-CUETP8M1-38T-PU1p3")==0) dir_input+="PYTHIA8-CUETP8M1-38T/";
  if(strcmp(type,"HERWIGPP-CUETHS1-38T-PU1p3")==0) dir_input+="HERWIGPP-CUETHS1-38T/";
  if(strcmp(type,"EPOS-38T-PU1p3")==0)             dir_input+="EPOS-38T/";              
  
  MainAnalyzer* m = new MainAnalyzer();

  //-- track analysis
  if (strcmp(analysis,"track") == 0) {

    //-- output directory
    if(strcmp(selection,"nominal") == 0) dir_output+="-nominal/";			 
    if(strcmp(selection,"HFup") == 0)    dir_output+="-HFup/";			    	 
    if(strcmp(selection,"HFdown") == 0)  dir_output+="-HFdown/";                        
    
    if(strcmp(selection,"no_pu_reweight") == 0) dir_output+="-no-pu-reweight/";	    
    if(strcmp(selection,"no_pixel_cut") == 0)   dir_output+="-no-pixel-cut/";	    
    if(strcmp(selection,"no_vz_reweight") == 0) dir_output+="-no-vz-reweight/";	    

    if(strcmp(selection,"combined") == 0) dir_output+="-combined/";	    
    
    //-- data 13 TeV 38T	
    if (strcmp(type,"data-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data at 38T with "<<selection<<endl;      
      CreateDirectory(dir,dir_output);
      m->makeHistoTrack(dir_input,"data_Run2015B_ZeroBias",dir_output,type,selection);
    }
    
    //-- data 13 TeV HFRereco 38T	
    if (strcmp(type,"data-HFRereco-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data with HFRereco at 38T with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoTrack(dir_input,"data_Run2015B_ZeroBias",dir_output,type,selection);
    }
    
    //-- data 13 TeV PromptReco 38T	
    if (strcmp(type,"data-PromptReco-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data with PromptReco at 38T with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoTrack(dir_input,"data_Run2015B_ZeroBias",dir_output,type,selection);
    }
    
    //-- PYTHIA8-CUETP8M1 13 TeV 38T no PU                                                        
    if (strcmp(type,"PYTHIA8-CUETP8M1-38T-noPU")==0) {
      cout<<endl<<"We'll process 13 TeV PYTHIA8-CUETP8M1 at 38T no PU with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoTrack(dir_input,"PYTHIA8-CUETP8M1-38T-noPU",dir_output,type,selection);
    }
    
    //-- PYTHIA8-CUETP8M1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"PYTHIA8-CUETP8M1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV PYTHIA8-CUETP8M1 at 38T PU 1p3 with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoTrack(dir_input,"PYTHIA8-CUETP8M1-38T-PU1p3",dir_output,type,selection);
    }
    
    //-- HERWIGPP-CUETHS1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"HERWIGPP-CUETHS1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV HERWIGPP-CUETHS1 at 38T PU 1p3 with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoTrack(dir_input,"HERWIGPP-CUETHS1-38T-PU1p3",dir_output,type,selection);
    }
    
    //-- EPOS 13 TeV 38T PU 1p3                                                                                                                
    if (strcmp(type,"EPOS-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV EPOS at 38T PU 1p3 with "<<selection<<endl;
      CreateDirectory(dir,dir_output);      
      m->makeHistoTrack(dir_input,"EPOS-38T-PU1p3",dir_output,type,selection);
    }

  } //-- end track analysis

  //-- correlation analysis
  if (strcmp(analysis,"correlation") == 0) {

    //-- output directory
    dir_output+="-rho";

    if(strcmp(selection,"nominal") == 0) dir_output+="-nominal/";			 
    if(strcmp(selection,"HFup") == 0)    dir_output+="-HFup/";			    	 
    if(strcmp(selection,"HFdown") == 0)  dir_output+="-HFdown/";                        
    
    if(strcmp(selection,"no_pu_reweight") == 0) dir_output+="-no-pu-reweight/";	    
    if(strcmp(selection,"no_pixel_cut") == 0)   dir_output+="-no-pixel-cut/";	    
    if(strcmp(selection,"no_vz_reweight") == 0) dir_output+="-no-vz-reweight/";	    
    
    //-- data 13 TeV 38T	
    if (strcmp(type,"data-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data at 38T with "<<selection<<endl;      
      CreateDirectory(dir,dir_output);
      m->makeHistoRho(dir_input,"data_Run2015B_ZeroBias",dir_output,type,selection);
    }
    
    //-- data 13 TeV HFRereco 38T	
    if (strcmp(type,"data-HFRereco-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data with HFRereco at 38T with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoRho(dir_input,"data_Run2015B_ZeroBias",dir_output,type,selection);
    }
    
    //-- data 13 TeV PromptReco 38T	
    if (strcmp(type,"data-PromptReco-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data with PromptReco at 38T with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoRho(dir_input,"data_Run2015B_ZeroBias",dir_output,type,selection);
    }
    
    //-- PYTHIA8-CUETP8M1 13 TeV 38T no PU                                                        
    if (strcmp(type,"PYTHIA8-CUETP8M1-38T-noPU")==0) {
      cout<<endl<<"We'll process 13 TeV PYTHIA8-CUETP8M1 at 38T no PU with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoRho(dir_input,"PYTHIA8-CUETP8M1-38T-noPU",dir_output,type,selection);
    }
    
    //-- PYTHIA8-CUETP8M1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"PYTHIA8-CUETP8M1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV PYTHIA8-CUETP8M1 at 38T PU 1p3 with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoRho(dir_input,"PYTHIA8-CUETP8M1-38T-PU1p3",dir_output,type,selection);
    }
    
    //-- HERWIGPP-CUETHS1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"HERWIGPP-CUETHS1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV HERWIGPP-CUETHS1 at 38T PU 1p3 with "<<selection<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoRho(dir_input,"HERWIGPP-CUETHS1-38T-PU1p3",dir_output,type,selection);
    }
    
    //-- EPOS 13 TeV 38T PU 1p3                                                                                                                
    if (strcmp(type,"EPOS-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV EPOS at 38T PU 1p3 with "<<selection<<endl;
      CreateDirectory(dir,dir_output);      
      m->makeHistoRho(dir_input,"EPOS-38T-PU1p3",dir_output,type,selection);
    }

  } //-- end correlation analysis

  //-- vertex analysis
  if (strcmp(analysis,"vertex") == 0) {

    //-- output directory
    dir_output+="-vertex-weight/";			 
   												    
    //-- data 13 TeV 38T	
    if (strcmp(type,"data-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data at 38T with iteration "<<iteration<<endl;      
      CreateDirectory(dir,dir_output);
      m->makeHistoVertex(dir_input,"data_Run2015B_ZeroBias",dir_output,type,iteration);
    }
    
    //-- PYTHIA8-CUETP8M1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"PYTHIA8-CUETP8M1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV PYTHIA8-CUETP8M1 at 38T PU 1p3 with iteration "<<iteration<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoVertex(dir_input,"PYTHIA8-CUETP8M1-38T-PU1p3",dir_output,type,iteration);
    }
    
    //-- HERWIGPP-CUETHS1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"HERWIGPP-CUETHS1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV HERWIGPP-CUETHS1 at 38T PU 1p3 with iteration "<<iteration<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoVertex(dir_input,"HERWIGPP-CUETHS1-38T-PU1p3",dir_output,type,iteration);
    }
    
    //-- EPOS 13 TeV 38T PU 1p3                                                                                                                
    if (strcmp(type,"EPOS-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV EPOS at 38T PU 1p3 with iteration "<<iteration<<endl;
      CreateDirectory(dir,dir_output);      
      m->makeHistoVertex(dir_input,"EPOS-38T-PU1p3",dir_output,type,iteration);
    }
  } //-- end vertex analysis

  //-- pu analysis
  if (strcmp(analysis,"pu") == 0) {

    //-- output directory
    dir_output+="-pu-weight/";			 
   												    
    //-- data 13 TeV 38T	
    if (strcmp(type,"data-38T") == 0) {
      cout<<endl<<"We'll process 13 TeV data at 38T with iteration "<<iteration<<endl;      
      CreateDirectory(dir,dir_output);
      m->makeHistoPU(dir_input,"data_Run2015B_ZeroBias",dir_output,type,iteration);
    }
    
    //-- PYTHIA8-CUETP8M1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"PYTHIA8-CUETP8M1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV PYTHIA8-CUETP8M1 at 38T PU 1p3 with iteration "<<iteration<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoPU(dir_input,"PYTHIA8-CUETP8M1-38T-PU1p3",dir_output,type,iteration);
    }
    
    //-- HERWIGPP-CUETHS1 13 TeV 38T PU 1p3                                                                                                              
    if (strcmp(type,"HERWIGPP-CUETHS1-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV HERWIGPP-CUETHS1 at 38T PU 1p3 with iteration "<<iteration<<endl;
      CreateDirectory(dir,dir_output);
      m->makeHistoPU(dir_input,"HERWIGPP-CUETHS1-38T-PU1p3",dir_output,type,iteration);
    }
    
    //-- EPOS 13 TeV 38T PU 1p3                                                                                                                
    if (strcmp(type,"EPOS-38T-PU1p3")==0) {
      cout<<endl<<"We'll process 13 TeV EPOS at 38T PU 1p3 with iteration "<<iteration<<endl;
      CreateDirectory(dir,dir_output);      
      m->makeHistoPU(dir_input,"EPOS-38T-PU1p3",dir_output,type,iteration);
    }
  } //-- end pu analysis
 
  delete m;
    
  return(0);
}


