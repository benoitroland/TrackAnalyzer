#ifndef TrackAnalyzer_h
#define TrackAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "TrackCommon.h"

class TrackAnalyzer: public TrackCommon {
  
 public:
  
  TrackAnalyzer();
  virtual ~TrackAnalyzer();
  
  void Loop(TString inputdir, TObjArray* filelist, TString type, TString selection, TString file_name);

 private:
};

#endif
