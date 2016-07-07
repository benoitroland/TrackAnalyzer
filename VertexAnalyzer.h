#ifndef VertexAnalyzer_h
#define VertexAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "VertexCommon.h"

class VertexAnalyzer: public VertexCommon {
  
 public:
  
  VertexAnalyzer();
  virtual ~VertexAnalyzer();
  
  void Loop(TString inputdir, TObjArray* filelist, TString type, int iteration, TString file_name);

 private:
};

#endif
