#ifndef MainAnalyzer_h
#define MainAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TCanvas.h>

#include "FileReader.h"

#include "TrackAnalyzer.h"
#include "RhoAnalyzer.h"
#include "VertexAnalyzer.h"
#include "PUAnalyzer.h"

class MainAnalyzer {
	public:
        MainAnalyzer();
		virtual ~MainAnalyzer();

		void makeHistoTrack(TString inputdir, TString name, TString outputdir, TString type, TString selection);
		void makeHistoRho(TString inputdir, TString name, TString outputdir, TString type, TString selection);
		void makeHistoVertex(TString inputdir, TString name, TString outputdir, TString type, int iteration);
		void makeHistoPU(TString inputdir, TString name, TString outputdir, TString type, int iteration);

		void saveAllCanvas(TString inputdir, TString name);
		void saveAllCanvasPDF(TString inputdir,TString name);

		void setPlotStyle();
		void setCMSStyle(); 

	private:
		FileReader reader_;

		TrackAnalyzer trackAnalyzer_;
		RhoAnalyzer rhoAnalyzer_;
		VertexAnalyzer vertexAnalyzer_;
		PUAnalyzer puAnalyzer_;

		std::vector<TCanvas*> canvasvector_;
};

#endif
