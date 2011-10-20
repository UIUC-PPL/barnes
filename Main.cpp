#include "Main.h"
#include "Parameters.h"

#include "TreePiece.h"
#include "DataManager.h"
#include "Reduction.h"
#include "Messages.h"

#include "defaults.h"

#include <iostream>
#include <fstream>
using namespace std;

CProxy_TreePiece treePieceProxy;
CProxy_DataManager dataManagerProxy;
CProxy_Main mainProxy;

Parameters globalParams;

Main::Main(CkArgMsg *msg){
  setParameters(msg);

  dataManagerProxy = CProxy_DataManager::ckNew();

  CkArrayOptions opts(globalParams.numTreePieces);
  /*
  CProxy_RRMap myMap = CProxy_RRMap::ckNew();
  opts.setMap(myMap);
  */
  treePieceProxy = CProxy_TreePiece::ckNew(opts);

#ifdef TRACE_REMOTE_DATA_REQUESTS
  traceRegisterUserEvent("Node",REMOTE_NODE_REQUEST);
  traceRegisterUserEvent("Particles",REMOTE_PARTICLE_REQUEST);
#endif

  mainProxy = thisProxy;

  thisProxy.commence();

  CkCallback cb(CkIndex_Main::quiescence(),thisProxy);
  CkStartQD(cb);

  numQuiescenceRecvd = 0;
  delete msg;
}

void Main::setParameters(CkArgMsg *m){
  map<string,string> table;
  params.extractParameters(m->argc, m->argv, table); 

  map<string,string>::iterator it = table.find("in");
  if(it == table.end()){
    CkPrintf("[Main] Must specify input particle file:\n");
    usage();
    CkAbort("bad command line arguments\n");
  }

  globalParams.filename = params.getsparam("in", table);

  globalParams.theta = params.getrparam("theta", DEFAULT_THETA, table);
  CkPrintf("theta: %f\n", globalParams.theta);

  globalParams.dtime = params.getrparam("dtime", DEFAULT_DTIME, table);
  CkPrintf("dtime: %f\n", globalParams.dtime);

  globalParams.dthf = globalParams.dtime/2.0;

  globalParams.epssq = params.getrparam("eps", DEFAULT_EPS, table);
  CkPrintf("eps: %f\n", globalParams.epssq);
  globalParams.epssq = globalParams.epssq*globalParams.epssq;

  globalParams.tolsq = globalParams.theta;
  CkPrintf("tol: %f\n", globalParams.tolsq);
  globalParams.tolsq = globalParams.tolsq*globalParams.tolsq;

  globalParams.ppc = params.getiparam("ppc", DEFAULT_PPC, table); 
  CkPrintf("ppc: %d\n", globalParams.ppc);

  globalParams.ppb = params.getiparam("b", DEFAULT_PPB, table);
  CkPrintf("bucketSize: %d\n", globalParams.ppb);

  globalParams.iterations = params.getiparam("killat", DEFAULT_KILLAT, table);
  CkPrintf("killat: %d\n", globalParams.iterations);
  
  globalParams.cacheLineSize = params.getiparam("chunkdepth", DEFAULT_CHUNK_DEPTH, table);
  CkPrintf("chunkdepth: %d\n", globalParams.cacheLineSize);

  globalParams.yieldPeriod = params.getiparam("yield", DEFAULT_YIELD_PERIOD, table);
  CkPrintf("yieldPeriod: %d\n", globalParams.yieldPeriod);

  globalParams.balancePeriod = params.getiparam("balancePeriod", DEFAULT_BALANCE_PERIOD, table);
  CkPrintf("balancePeriod: %d\n", globalParams.balancePeriod);

  getNumParticles();

  it = table.find("p");
  if(it == table.end()){
    globalParams.numTreePieces = ((Real)globalParams.numParticles/((Real)globalParams.ppc))*2.0;
    if(globalParams.numTreePieces == 0) globalParams.numTreePieces = 1;
  }
  else{
    globalParams.numTreePieces = atoi(it->second.c_str());
  }

  CkPrintf("tree pieces: %d\n", globalParams.numTreePieces);


}

void Main::getNumParticles(){
  ifstream partFile;
  CkPrintf("[Main] file %s\n", globalParams.filename.c_str());
  partFile.open(globalParams.filename.c_str(), ios::in | ios::binary);
  CkAssert(partFile.is_open());

  partFile.read((char *)(&globalParams.numParticles),sizeof(int)); 
  CkPrintf("[Main] numParticles %d\n", globalParams.numParticles);
  partFile.close();
}

void Main::commence(){
  CkPrintf("[Main] load particles\n");
  CkReductionMsg *redMsg;
  dataManagerProxy.loadParticles(CkCallbackResumeThread((void *&)redMsg));

  BoundingBox &universe = *((BoundingBox *)redMsg->getData());
  Real pad = 0.00001;
  universe.expand(pad);
  dataManagerProxy.decompose(universe);

}

void Main::niceExit(){
  CkPrintf("[Main] graceful exit\n");
  CkExit();
}

void Main::quiescence(){
  CkPrintf("[Main] Quiescence detected!\n");
  treePieceProxy.quiescence();
  dataManagerProxy.quiescence();
}

void Main::quiescenceExit(){
  numQuiescenceRecvd++;
  if(numQuiescenceRecvd == 2){
    CkExit();
  }
}

string NodeTypeString[] = { 
  "Invalid",
  "Internal",
  "Bucket",
  "EmptyBucket",
  "Boundary",
  "Remote",
  "RemoteBucket",
  "RemoteEmptyBucket"
};

string NodeTypeColor[] = {
  "firebrick1",
  "darkolivegreen1",
  "darkolivegreen3",
  "darksalmon",
  "darkkhaki",
  "deepskyblue1",
  "dodgerblue4",
  "deeppink"
};

#if 0 
#endif

void Main::usage(){
  map<string,string> usage;
  usage["in"] = "input file";
  usage["ppc"] = "particles per chare";
  usage["b"] = "particles per bucket (leaf)";
  usage["theta"] = "opening angle";
  usage["killat"] = "num single steps";
  usage["chunkDepth"] = "when fetching remote data, what depth of subtree to fetch";
  usage["yield"] = "how many buckets to process before yielding processor";
  usage["ppc"] = "particleschare";


  map<string,string>::iterator it;
  CkPrintf("Usage: ./charmrun +p<numproc> ./barnes <options>\n");
  CkPrintf("Options take the form -<name>=<value>\n");
  CkPrintf("Options are described below:\n");
  for(it = usage.begin(); it != usage.end(); it++){
    CkPrintf("%s : %s\n", it->first.c_str(), it->second.c_str());
  }
}

#include "barnes.def.h"
