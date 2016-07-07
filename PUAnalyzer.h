#ifndef PUAnalyzer_h
#define PUAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "PUCommon.h"

class PUAnalyzer: public PUCommon {
  
 public:
  
  PUAnalyzer();
  virtual ~PUAnalyzer();
  
  void Loop(TString inputdir, TObjArray* filelist, TString type, int iteration, TString file_name);

 private:
};

#endif
