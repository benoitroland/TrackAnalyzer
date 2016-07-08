#ifndef Common_h
#define Common_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>

class Common {

 public:

  Common();
  virtual ~Common();

  void CheckHisto(TH1D* h);
  void CheckHisto(TH2D* h);

  void CheckGraph(TGraphErrors* graph);

  void SetAxisName(TH1D* h,TString xleg, TString yleg);

  TH1D* MakeHisto(TString name, TString title, TString xleg, TString yleg, int nbin, double bmin, double bmax);
  TH1D* MakeHisto(TString name, TString title, TString xleg, TString yleg, int nbin, double* listbin);

  TH2D* MakeHisto(TString name, TString title, TString xleg, TString yleg, int nbinx, double bxmin, double bxmax, int nbiny, double bymin, double bymax);
  TH2D* MakeHisto(TString name, TString title, TString xleg, TString yleg, int nbinx, double* listbinx, int nbiny, double* listbiny);

  void MakeGraph(TGraphErrors* graph, TString name, TString title, TString xleg, TString yleg);
  TProfile* MakeProfile(TString name, TString title, TString xleg, TString yleg, int nbin, double bmin, double bmax);

  void ComputeRatio(double a, double b, double error_a, double error_b, double &ratio, double &error);
  void ComputeProduct(double a, double b, double error_a, double error_b, double &product, double &error);

  void ComputeHistoRatio(TH1D* hup,TH1D* hdown, TH1D& hratio);
  void ComputeHistoProduct(TH1D* ha, TH1D* hb, TH1D& hproduct);

  void ComputeOneMinusFake(TH1D* hFake, TH1D& hOneMinusFake);
  void ComputeCFPt(TH1D* hOneMinusBkg, TH1D* hAcc, TH1D* hPurity, TH1D* hStability, TH1D& hCF);

  int GetPhiMapping(int sector);

  void GetHistoContent(TH1D* histo, double &content, double &error);

  double GetDeltaR(double phi_trackjet, double eta_trackjet, double phi_genjet, double eta_genjet);

  int GetSectorSubtracted(int isector_max, int shift);
  int GetSectorAdded(int isector_max, int shift);

  double GetDeltaPhi(double phi1, double phi2);
  double Rad2Deg(double phi);

  TF1* FitGauss(TH1D* histo,double min, double max, double DeltaM);

  bool IsGoodTrack(int quality,double pt, double pterror, double dxy, double dxyerror, double dz, double dzerror, int npixelhit);

  double GetVertexWeight(double Vz,bool* AllMC, int iteration);
  void GetPUWeight(bool* AllMC,double* weight,int iteration);

  bool isProton(bool* AllMC, int id);
  bool isPion(bool* AllMC, int id);
  bool isKaon(bool* AllMC, int id);

 private:
};

#endif

