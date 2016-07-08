#ifndef RhoAnalyzer_h
#define RhoAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "RhoCommon.h"

class RhoAnalyzer: public RhoCommon {
  
 public:
  
  RhoAnalyzer();
  virtual ~RhoAnalyzer();
  
  void Loop(TString inputdir, TObjArray* filelist, TString type, TString selection, TString file_name);
  

 private:
};

#endif
