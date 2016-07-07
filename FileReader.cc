#include "FileReader.h"

#include "TObjString.h"
#include "TSystem.h"
#include <iostream>

FileReader::FileReader() {

}

FileReader::~FileReader() { 

}

TObjArray* FileReader::getFileList(TString inputdir, TString regexpstr) {

	TObjArray*  filelist = new TObjArray;
	filelist->Clear();
	void *dir = gSystem->OpenDirectory(inputdir);
	if (!dir) {
		std::cout <<  "Error in FileReader: couldn't open directory" << std::endl;
		return 0;
	}
		
	const char* file;
	TString     fileName;
	TRegexp     regexp(regexpstr);
	while ((file = gSystem->GetDirEntry(dir))) {
		fileName = file;
		if (!fileName.CompareTo(".") || !fileName.CompareTo(".."))
			continue;
		if (fileName.Index(regexp) != kNPOS)
			filelist->Add(new TObjString(fileName));
	}  
	gSystem->FreeDirectory(dir);
	
	TIter       next(filelist); 
	TObjString* fn = 0;
	
	// do some print out
	while((fn = (TObjString*)next())) { 
	  std::cout << "FileReader found file: " << fn->GetString() << " matching your criteria"<<std::endl; // in directory " << inputdir << std::endl;
	}
	
	return filelist;

}
