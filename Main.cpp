/*
  CharmBH - Main.cpp

  Main is the mainchare; control starts here. Main initializes the computation
  by reading input parameters and initiates particle loading.
  Once they have been loaded, the mainchare directs all PEs to begin
  the Oct-decompostion of particles onto PEs.
*/

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

/* Proxies for tree piece array and data manager group */
CProxy_TreePiece treePieceProxy;
CProxy_DataManager dataManagerProxy;
CProxy_Main mainProxy;

Parameters globalParams;

Main::Main(CkArgMsg *msg){
  /* Read input parameters from command line. Fill in defaults for
     unspecified parameters.
  */
  setParameters(msg);

  /* Create the DataManager group; there is one representative/member
     of the group on each PE.
  */
  dataManagerProxy = CProxy_DataManager::ckNew();

  CkArrayOptions opts(globalParams.numTreePieces);
  /* 
    We would like the tree piece array elements to be mapped to PEs
    in a round-robin fashion to begin with.
  */
  CProxy_RRMap myMap = CProxy_RRMap::ckNew();
  opts.setMap(myMap);
  /*
    Create the tree piece chare array.
  */
  treePieceProxy = CProxy_TreePiece::ckNew(opts);

#ifdef TRACE_REMOTE_DATA_REQUESTS
  traceRegisterUserEvent("Node",REMOTE_NODE_REQUEST);
  traceRegisterUserEvent("Particles",REMOTE_PARTICLE_REQUEST);
#endif

  /*
    Set mainProxy readonly so that it is available
    on all other PEs after this method finishes.
  */
  mainProxy = thisProxy;

  /*
    Begin actual computation
  */
  thisProxy.commence();

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

  /* Input file */
  string filename = params.getsparam("in", table);
  globalParams.nchars = filename.length();
  globalParams.filename = new char[globalParams.nchars+1];
  filename.copy(globalParams.filename,globalParams.nchars);
  globalParams.filename[globalParams.nchars] = '\0';

  /* Opening angle */
  globalParams.theta = params.getrparam("theta", DEFAULT_THETA, table);
  CkPrintf("theta: %f\n", globalParams.theta);

  /* Step time */
  globalParams.dtime = params.getrparam("dtime", DEFAULT_DTIME, table);
  CkPrintf("dtime: %f\n", globalParams.dtime);

  /* Half of step time; used in leapfrog "kick" */
  globalParams.dthf = globalParams.dtime/2.0;

  /* Softening parameter */
  globalParams.epssq = params.getrparam("eps", DEFAULT_EPS, table);
  CkPrintf("eps: %f\n", globalParams.epssq);
  globalParams.epssq = globalParams.epssq*globalParams.epssq;

  globalParams.tolsq = globalParams.theta;
  CkPrintf("tol: %f\n", globalParams.tolsq);
  globalParams.tolsq = globalParams.tolsq*globalParams.tolsq;

  /* Particles per chare (number of particles per tree piece) */
  globalParams.ppc = params.getiparam("ppc", DEFAULT_PPC, table); 
  CkPrintf("ppc: %d\n", globalParams.ppc);

  globalParams.ppb = params.getiparam("b", DEFAULT_PPB, table);
  CkPrintf("bucketSize: %d\n", globalParams.ppb);

  /* Number of steps to run the simulation for */
  globalParams.iterations = params.getiparam("killat", DEFAULT_KILLAT, table);
  CkPrintf("killat: %d\n", globalParams.iterations);
  
  /* 
    When a remote node is requested, the subtree beneath
    it up to "cacheLineSize" is communicated alongwith the
    node itself. This helps to coarsen the grain of communication.
  */
  globalParams.cacheLineSize = params.getiparam("chunkdepth", DEFAULT_CHUNK_DEPTH, table);
  CkPrintf("chunkdepth: %d\n", globalParams.cacheLineSize);

  /*
    When traversing the global tree for each bucket of particles, 
    do "yield" number of buckets before returning to the scheduler
    to look for outstanding request messages, etc.
  */
  globalParams.yieldPeriod = params.getiparam("yield", DEFAULT_YIELD_PERIOD, table);
  CkPrintf("yieldPeriod: %d\n", globalParams.yieldPeriod);

  /*
    Do load balancing every "balancePeriod" iterations.
  */
  globalParams.balancePeriod = params.getiparam("balancePeriod", DEFAULT_BALANCE_PERIOD, table);
  CkPrintf("balancePeriod: %d\n", globalParams.balancePeriod);

  globalParams.doPrintAccel = false;
  if(table.find("output") != table.end()) globalParams.doPrintAccel = true;
  CkPrintf("doPrintAccel: %d\n", globalParams.doPrintAccel);

  /* Number of steps between decompositions */
  globalParams.decompPeriod = params.getiparam("decompPeriod", DEFAULT_DECOMP_PERIOD, table); 
  CkPrintf("decompPeriod: %d\n", globalParams.decompPeriod);

  /* To prevent sending of very large messages after decomposition */
  globalParams.particleMsgMaxSize = params.getiparam("particleMsgMaxSize", DEFAULT_PARTICLE_MSG_MAX_SIZE, table); 
  CkPrintf("particleMsgMaxSize: %d\n", globalParams.particleMsgMaxSize);
  CkAssert(globalParams.particleMsgMaxSize > sizeof(Particle));


  /*
    Use the filename obtained previously to read just the 
    number of particles from the input file.
  */
  getNumParticles();

  /* 
    Find the appropriate number of tree pieces to generate.
  */
  it = table.find("p");
  /* 
    Find the appropriate number of tree pieces to generate.
  */
  if(it == table.end()){
    /*
      If the number of pieces is not specified by the user, 
      calculate it from the number of particles per chare and
      the total number of particles.
    */
    globalParams.numTreePieces = ((Real)globalParams.numParticles/((Real)globalParams.ppc))*2.0;
    if(globalParams.numTreePieces == 0) globalParams.numTreePieces = 1;
  }
  else{
    /*
      If the number of pieces is specified by the user, use that.
    */
    globalParams.numTreePieces = atoi(it->second.c_str());
  }

  CkPrintf("tree pieces: %d\n", globalParams.numTreePieces);


}

void Main::getNumParticles(){
  ifstream partFile;
  CkPrintf("[Main] file %s\n", globalParams.filename);
  partFile.open(globalParams.filename, ios::in | ios::binary);
  CkAssert(partFile.is_open());

  partFile.read((char *)(&globalParams.numParticles),sizeof(int)); 
  CkPrintf("[Main] numParticles %d\n", globalParams.numParticles);
  partFile.close();
}

void Main::commence(){
  CkPrintf("[Main] load particles\n");
  CkReductionMsg *redMsg;
  /* 
    Tell each PE to read the particles 
    from its portion of the input file
  */
  dataManagerProxy.loadParticles(CkCallbackResumeThread((void *&)redMsg));

  /*
    Each PE contributes the bounding box of the particles
    it read to a reduction; as a result, the bounding box
    of all particles in the simulated "universe" is obtained.
  */
  BoundingBox &universe = *((BoundingBox *)redMsg->getData());
  Real pad = 0.00001;
  /*
    Provide a little padding so that the particles closest to
    the extremities of the universe get correct keys.
  */
  universe.expand(pad);
  /*
    Tell all PEs to begin Oct decomposition of read particles.
  */
  
  CkCallback cb(CkIndex_DataManager::processSubmittedParticles(),dataManagerProxy);
  CkStartQD(cb);
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
