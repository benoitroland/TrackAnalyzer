#ifndef FileReader_h
#define FileReader_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>

class FileReader {
	public:
        FileReader();
	virtual ~FileReader();
        TObjArray* getFileList(TString inputdir, TString regexpstr);

	private:
};

#endif
