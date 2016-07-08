{
  gROOT->ProcessLine(".L ./FileReader.cc+");

  gROOT->ProcessLine(".L ./Common.cc+");

  gROOT->ProcessLine(".L ./TrackCommon.cc+");
  gROOT->ProcessLine(".L ./TrackAnalyzer.cc+");

  gROOT->ProcessLine(".L ./VertexCommon.cc+");  
  gROOT->ProcessLine(".L ./VertexAnalyzer.cc+");

  gROOT->ProcessLine(".L ./PUCommon.cc+");  
  gROOT->ProcessLine(".L ./PUAnalyzer.cc+");

  gROOT->ProcessLine(".L ./RhoCommon.cc+");
  gROOT->ProcessLine(".L ./RhoAnalyzer.cc+");
  
  gROOT->ProcessLine(".L ./MainAnalyzer.cc+");
}
