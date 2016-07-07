SHELL = /bin/bash
COMPILER = g++ -g
COMPILERSO = g++ -g +a1 -b
LIBDIR = ${PWD}
CFLAGS = $(shell root-config --cflags)
LIBS = $(shell root-config --libs) -lPhysics -lThread -lMinuit -lHtml -lVMC -lEG -lGeom -Wl,-rpath -Wl,${LIBDIR}
ROOT = /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_2/external/slc6_amd64_gcc491/bin/root -b -q

COMMON = Common.o Common_cc.so
TRACKANA = TrackCommon.o TrackCommon_cc.so TrackAnalyzer.o TrackAnalyzer_cc.so
RHOANA = RhoCommon.o RhoCommon_cc.so RhoAnalyzer.o RhoAnalyzer_cc.so
VERTEXANA = VertexCommon.o VertexCommon_cc.so VertexAnalyzer.o VertexAnalyzer_cc.so
PUANA = PUCommon.o PUCommon_cc.so PUAnalyzer.o PUAnalyzer_cc.so

.SUFFIXES: .cc .o

all : Run 

Run : Run.cc MainAnalyzer.o MainAnalyzer_cc.so ${COMMON} ${TRACKANA} ${RHOANA} ${VERTEXANA} ${PUANA} FileReader.o FileReader_cc.so 
	${COMPILER} Run.cc -o $@ MainAnalyzer.o MainAnalyzer_cc.so ${COMMON} ${TRACKANA} ${RHOANA} ${VERTEXANA} ${PUANA} FileReader.o FileReader_cc.so ${CFLAGS} ${LIBS}


MainAnalyzer_cc.so :
	${ROOT} Build.cc+

MainAnalyzer.o : MainAnalyzer.cc MainAnalyzer.h FileReader.cc FileReader.h 
	${COMPILER} -c MainAnalyzer.cc FileReader.cc  ${CFLAGS}

Common.o : Common.cc Common.h   
	${COMPILER} -c Common.cc ${CFLAGS}

TrackCommon.o : TrackCommon.cc TrackCommon.h   
	${COMPILER} -c TrackCommon.cc ${CFLAGS}

TrackAnalyzer.o : TrackAnalyzer.cc TrackAnalyzer.h
	${COMPILER} -c TrackAnalyzer.cc ${CFLAGS}  

RhoCommon.o : RhoCommon.cc RhoCommon.h   
	${COMPILER} -c RhoCommon.cc ${CFLAGS}

RhoAnalyzer.o : RhoAnalyzer.cc RhoAnalyzer.h
	${COMPILER} -c RhoAnalyzer.cc ${CFLAGS}  

VertexCommon.o : VertexCommon.cc VertexCommon.h
	${COMPILER} -c VertexCommon.cc ${CFLAGS}

VertexAnalyzer.o : VertexAnalyzer.cc VertexAnalyzer.h
	${COMPILER} -c VertexAnalyzer.cc ${CFLAGS}    

PUCommon.o : PUCommon.cc PUCommon.h	      
	${COMPILER} -c PUCommon.cc ${CFLAGS}      

PUAnalyzer.o : PUAnalyzer.cc PUAnalyzer.h 
	${COMPILER} -c PUAnalyzer.cc ${CFLAGS}    

FileReader.o : FileReader.cc FileReader.h
	${COMPILER} -c FileReader.cc ${CFLAGS}

clean :
	rm *.pcm *.o *.so *_cc.d Run

