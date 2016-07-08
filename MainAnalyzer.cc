#include "MainAnalyzer.h"

#include <TROOT.h>
#include "TObjString.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TList.h"
#include "TKey.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TSystem.h"
#include <TStyle.h>

#include <iostream>

using namespace std;

MainAnalyzer::MainAnalyzer() {
	
	// initialize helper classes
	//reader_ = new FileReader();
	//histogetter_ = new HistoRetriever();
	//treeanalyzer_ = new TreeAnalyzer();

}

MainAnalyzer::~MainAnalyzer() { 
	

}


void MainAnalyzer::makeHistoTrack(TString inputdir, TString name, TString outputdir, TString type, TString selection) {

  TObjArray *inputfile = reader_.getFileList(inputdir,name);
  
  TString outputfile = TString(outputdir) + TString("output_trackanalysis_") + TString(type) + TString("_") + TString(selection) + TString(".root");

  trackAnalyzer_.Loop(inputdir,inputfile,type,selection,outputfile);

  cout<<"output file name: "<<outputfile<<endl;

}

void MainAnalyzer::makeHistoRho(TString inputdir, TString name, TString outputdir, TString type, TString selection) {

  TObjArray *inputfile = reader_.getFileList(inputdir,name);
  
  TString outputfile = TString(outputdir) + TString("output_rhoanalysis_") + TString(type) + TString("_") + TString(selection) + TString(".root");

  rhoAnalyzer_.Loop(inputdir,inputfile,type,selection,outputfile);

  cout<<"output file name: "<<outputfile<<endl;

}

void MainAnalyzer::makeHistoVertex(TString inputdir, TString name, TString outputdir, TString type, int iteration) {

  TObjArray *inputfile = reader_.getFileList(inputdir,name);
  
  TString outputfile = TString(outputdir) + TString("output_vertexanalysis_") + TString(type) + TString("_iteration_") + TString::Format("%d",iteration) + TString(".root");

  vertexAnalyzer_.Loop(inputdir,inputfile,type,iteration,outputfile);

  cout<<"output file name: "<<outputfile<<endl;

}

void MainAnalyzer::makeHistoPU(TString inputdir, TString name, TString outputdir, TString type, int iteration) {

  TObjArray *inputfile = reader_.getFileList(inputdir,name);

  TString outputfile = TString(outputdir) + TString("output_puanalysis_") + TString(type) + TString("_iteration_") + TString::Format("%d",iteration) + TString(".root");

  puAnalyzer_.Loop(inputdir,inputfile,type,iteration,outputfile);

  cout<<"output file name: "<<outputfile<<endl;

}

void MainAnalyzer::saveAllCanvasPDF(TString inputdir,TString name) {
  
  cout<<endl<<"begin to save canvas"<<endl;
  
  TString file_pdf = TString(inputdir) + TString("plot_") + TString(name) + TString(".pdf");
  
  canvasvector_[0]->Print(TString(TString(file_pdf)+TString("[")).Data());

  for (unsigned int i=0;i<canvasvector_.size();i++) {
    TCanvas *c = canvasvector_[i];
    canvasvector_[i]->Print(file_pdf.Data());
    TString cname;
    cname.Append(inputdir);
    cname.Append(c->GetName());
    cname.Append("_");
    cname.Append(name);
    cname.Append(".png");
    c->SaveAs(cname);
  }
  
  canvasvector_[0]->Print(TString(TString(file_pdf)+TString("]")).Data());
  cout<<"canvas saved !"<<endl;
}

void MainAnalyzer::saveAllCanvas(TString inputdir,TString name) {
  
  cout<<"begin to save canvas"<<endl;
  
  TString file_pdf = TString(inputdir) + TString("plot_") + TString(name) + TString(".pdf");
  
  canvasvector_[0]->Print(TString(TString(file_pdf)+TString("[")).Data());
  
  for (unsigned int i=0;i<canvasvector_.size();i++) {
    TCanvas *c = canvasvector_[i];
    canvasvector_[i]->Print(file_pdf.Data());
    TString cname;
    cname.Append(inputdir);
    cname.Append(c->GetName());
    cname.Append("_");
    cname.Append(name);
    cname.Append(".png");
    c->SaveAs(cname);
  }
  
  canvasvector_[0]->Print(TString(TString(file_pdf)+TString("]")).Data());
  cout<<"canvas saved !"<<endl;
}

void MainAnalyzer::setPlotStyle() {

  //-- set plot styles
  gStyle->SetOptStat(111111);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetStatColor(kWhite);
}


void MainAnalyzer::setCMSStyle(){
	
	std::cout << "CMS Style Loaded" << std::endl;
	
	TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
	
	// For the canvas:
	tdrStyle->SetCanvasBorderMode(0);
	tdrStyle->SetCanvasColor(kWhite);
	// tdrStyle->SetCanvasDefH(200); //Height of canvas
	//tdrStyle->SetCanvasDefW(200); //Width of canvas
	tdrStyle->SetCanvasDefX(0);   //POsition on screen
	tdrStyle->SetCanvasDefY(0);
	
	// For the Pad:
	tdrStyle->SetPadBorderMode(0);
	tdrStyle->SetPadBorderSize(0);
	tdrStyle->SetPadColor(kWhite);
	//   tdrStyle->SetPadGridX(false);
	//   tdrStyle->SetPadGridY(false);
	//   tdrStyle->SetGridColor(0);
	//   tdrStyle->SetGridStyle(3);
	//   tdrStyle->SetGridWidth(1);
	
	// For the frame:
	tdrStyle->SetFrameBorderMode(0);
	tdrStyle->SetFrameBorderSize(0);
	tdrStyle->SetFrameFillColor(0);
	tdrStyle->SetFrameFillStyle(0);
	tdrStyle->SetFrameLineColor(1);
	tdrStyle->SetFrameLineStyle(1);
	tdrStyle->SetFrameLineWidth(1);
	
	// For the histo:
	// tdrStyle->SetHistFillColor(1);
	// tdrStyle->SetHistFillStyle(0);
	tdrStyle->SetHistLineColor(1);
	tdrStyle->SetHistLineStyle(0);
	tdrStyle->SetHistLineWidth(1);
	// tdrStyle->SetLegoInnerR(Float_t rad =3D 0.5);
	// tdrStyle->SetNumberContours(Int_t number =3D 20);
	
	tdrStyle->SetEndErrorSize(2);
	//tdrStyle->SetErrorMarker(20);
	//tdrStyle->SetErrorX(0.);
	
	tdrStyle->SetMarkerStyle(20);
	
	//For the fit/function:
	tdrStyle->SetOptFit(0);
	tdrStyle->SetFitFormat("5.4g");
	tdrStyle->SetFuncColor(2);
	tdrStyle->SetFuncStyle(1);
	tdrStyle->SetFuncWidth(1);
	
	//For the date:
	tdrStyle->SetOptDate(0);
	// tdrStyle->SetDateX(Float_t x =3D 0.01);
	// tdrStyle->SetDateY(Float_t y =3D 0.01);
	
	// For the statistics box:
	//  tdrStyle->SetOptFile(0);
	tdrStyle->SetOptStat(0); // To display the mean and RMS:   = SetOptStat("mr");
	tdrStyle->SetStatColor(kWhite);
	tdrStyle->SetStatFont(42);
	tdrStyle->SetStatFontSize(0.025);
	tdrStyle->SetStatTextColor(kBlack);
	tdrStyle->SetStatFormat("6.4g");
	tdrStyle->SetStatBorderSize(0);
	tdrStyle->SetStatH(0.1);
	tdrStyle->SetStatW(0.15);
	// tdrStyle->SetStatStyle(Style_t style =3D 1001);
	// tdrStyle->SetStatX(Float_t x =3D 0);
	// tdrStyle->SetStatY(Float_t y =3D 0);
	
	// Margins:
	tdrStyle->SetPadTopMargin(0.05);
	tdrStyle->SetPadBottomMargin(0.13);
	tdrStyle->SetPadLeftMargin(0.13);
	tdrStyle->SetPadRightMargin(0.05);
	
	// For the Global title:
	
	tdrStyle->SetOptTitle(0);
	tdrStyle->SetTitleFont(42);
	tdrStyle->SetTitleColor(1);
	tdrStyle->SetTitleTextColor(kWhite);
	tdrStyle->SetTitleFillColor(kWhite);
	tdrStyle->SetTitleFontSize(0.05);
	// tdrStyle->SetTitleH(0); // Set the height of the title box
	// tdrStyle->SetTitleW(0); // Set the width of the title box
	// tdrStyle->SetTitleX(0); // Set the position of the title box
	// tdrStyle->SetTitleY(0.985); // Set the position of the title box
	// tdrStyle->SetTitleStyle(Style_t style =3D 1001);
	// tdrStyle->SetTitleBorderSize(2);
	
	// For the axis titles:
	
	tdrStyle->SetTitleColor(1, "XYZ");
	tdrStyle->SetTitleFont(42, "XYZ");
	tdrStyle->SetTitleSize(0.06, "XYZ");
	// tdrStyle->SetTitleXSize(Float_t size =3D 0.02); // Another way to = set the size?
	// tdrStyle->SetTitleYSize(Float_t size =3D 0.02);
	tdrStyle->SetTitleXOffset(0.9);
	tdrStyle->SetTitleYOffset(1.05);
	// tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the = Offset
	
	// For the axis labels:
	
	tdrStyle->SetLabelColor(1, "XYZ");
	tdrStyle->SetLabelFont(42, "XYZ");
	tdrStyle->SetLabelOffset(0.007, "XYZ");
	tdrStyle->SetLabelSize(0.05, "XYZ");
	
	// For the axis:
	
	tdrStyle->SetAxisColor(1, "XYZ");
	tdrStyle->SetStripDecimals(kTRUE);
	tdrStyle->SetTickLength(0.03, "XYZ");
	tdrStyle->SetNdivisions(510, "XYZ");
	tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side = of the frame
	tdrStyle->SetPadTickY(1);
	
	// Change for log plots:
	//   tdrStyle->SetOptLogx(0);
	//   tdrStyle->SetOptLogy(0);
	//   tdrStyle->SetOptLogz(0);
	
	// Postscript options:
	//  tdrStyle->SetPaperSize(7.5,7.5);
	// tdrStyle->SetLineScalePS(Float_t scale =3D 3);
	// tdrStyle->SetLineStyleString(Int_t i, const char* text);
	// tdrStyle->SetHeaderPS(const char* header);
	// tdrStyle->SetTitlePS(const char* pstitle);
	
	// tdrStyle->SetBarOffset(Float_t baroff =3D 0.5);
	// tdrStyle->SetBarWidth(Float_t barwidth =3D 0.5);
	// tdrStyle->SetPaintTextFormat(const char* format =3D "g");
	// tdrStyle->SetPalette(Int_t ncolors =3D 0, Int_t* colors =3D 0);
	// tdrStyle->SetTimeOffset(Double_t toffset);
	// tdrStyle->SetHistMinimumZero(kTRUE);
	
	tdrStyle->cd();
	// gSystem->Load("libRooStats");
	// using namespace RooFit;
	
}

