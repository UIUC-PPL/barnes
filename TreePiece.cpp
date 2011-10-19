#include "TreePiece.h"
#include "Messages.h"
#include "DataManager.h"
#include "Parameters.h"

#include <fstream>

extern CProxy_Main mainProxy;
extern CProxy_DataManager dataManagerProxy;
extern Parameters globalParams;

TreePiece::TreePiece() {
  usesAtSync = CmiTrue;
  iteration = 0;
  init();
}

void TreePiece::init(){
  myNumParticles = 0;
  numDecompMsgsRecvd = 0;
  decompMsgsRecvd.length() = 0;
  largestKey = Key(0);
  smallestKey = ~largestKey;
  totalNumTraversals = 2;
  myDM = dataManagerProxy.ckLocalBranch();
  myRoot = NULL;
  root = NULL;
  myBuckets = NULL;
  myNumBuckets = 0;
}

void TreePiece::receiveParticles(ParticleMsg *msg){
  int msgNumParticles = msg->numParticles;
  decompMsgsRecvd.push_back(msg);
  myNumParticles += msgNumParticles;
  numDecompMsgsRecvd++;
  if(smallestKey > msg->part[0].key) smallestKey = msg->part[0].key;
  if(largestKey < msg->part[msgNumParticles-1].key) largestKey = msg->part[msgNumParticles-1].key;
  
  if(numDecompMsgsRecvd == CkNumPes()){
    submitParticles();
    numDecompMsgsRecvd = 0;
  }
}

void TreePiece::receiveParticles(){
  numDecompMsgsRecvd++;
  if(numDecompMsgsRecvd == CkNumPes()){
    submitParticles();
    numDecompMsgsRecvd = 0;
  }
}

void TreePiece::submitParticles(){
  myDM->submitParticles(&decompMsgsRecvd,myNumParticles,this,smallestKey,largestKey);
}

void TreePiece::prepare(Node<ForceData> *_root, Node<ForceData> *_myRoot, Node<ForceData> **buckets, int bucketStart, int bucketEnd){
  root = _root;
  myRoot = _myRoot;
  myBuckets = buckets+bucketStart;
  myNumBuckets = bucketEnd-bucketStart;

  if(myRoot != NULL) CkPrintf("tree piece %d prepare root %lu pe %d \n", thisIndex, myRoot->getKey(), CkMyPe());
}

void TreePiece::startTraversal(){
  trav.setDataManager(myDM);
  
  if(myNumBuckets == 0){
    finishIteration();
    return;
  }

  remoteTraversalState.reset(this,myNumBuckets,myBuckets);
  remoteTraversalWorker.reset(this,&remoteTraversalState,*myBuckets);
  RescheduleMsg *msg = new (NUM_PRIORITY_BITS) RescheduleMsg;
  *(int *)CkPriorityPtr(msg) = REMOTE_GRAVITY_PRIORITY;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  thisProxy[thisIndex].doRemoteGravity(msg);

  localTraversalState.reset(this,myNumBuckets,myBuckets);
  localTraversalWorker.reset(this,&localTraversalState,*myBuckets);

  msg = new (NUM_PRIORITY_BITS) RescheduleMsg;
  *(int *)CkPriorityPtr(msg) = LOCAL_GRAVITY_PRIORITY;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  thisProxy[thisIndex].doLocalGravity(msg);
}

void TreePiece::doLocalGravity(RescheduleMsg *msg){
  int i;
  for(i = 0; i < globalParams.yieldPeriod && 
                 localTraversalState.current < myNumBuckets; 
                 i++){
    trav.topDownTraversal(root,&localTraversalWorker,&localTraversalState);
    localTraversalState.current++;
    localTraversalState.currentBucketPtr++;
    localTraversalWorker.setContext(*localTraversalState.currentBucketPtr);
  }

  localTraversalState.decrPending(i);
  if(localTraversalState.current < myNumBuckets){
    CkAssert(!localTraversalState.complete());
    thisProxy[thisIndex].doLocalGravity(msg);
  }
  else{
    delete msg;
    if(localTraversalState.complete()) traversalDone();
  }
}

void TreePiece::doRemoteGravity(RescheduleMsg *msg){
  int i;
  for(i = 0; i < globalParams.yieldPeriod &&
                 remoteTraversalState.current < myNumBuckets;
                 i++){
    trav.topDownTraversal(root,&remoteTraversalWorker,&remoteTraversalState);
    remoteTraversalState.current++;
    remoteTraversalState.currentBucketPtr++;
    remoteTraversalWorker.setContext(*remoteTraversalState.currentBucketPtr);
  }

  remoteTraversalState.decrPending(i);
  if(remoteTraversalState.current < myNumBuckets){
    CkAssert(!remoteTraversalState.complete());
    thisProxy[thisIndex].doRemoteGravity(msg);
  }
  else{
    delete msg;
    if(remoteTraversalState.complete()) traversalDone();
  }
}

void TreePiece::traversalDone(){ 
  totalNumTraversals--;
  if(totalNumTraversals == 0) finishIteration();
}

void TreePiece::finishIteration(){
  CmiUInt8 pn = localTraversalState.numInteractions[0]+remoteTraversalState.numInteractions[0];
  CmiUInt8 pp = localTraversalState.numInteractions[1]+remoteTraversalState.numInteractions[1];
  CmiUInt8 oc = localTraversalState.numInteractions[2]+remoteTraversalState.numInteractions[2];
  dataManagerProxy[CkMyPe()].traversalsDone(pn,pp,oc);

  localTraversalState.finishedIteration();
  remoteTraversalState.finishedIteration();
  
  init();

  iteration++;
}

void TreePiece::quiescence(){
  CkPrintf("QUIESCENCE tree piece %d proc %d submitted %d numBuckets %d trav_done %d outstanding local %d remote %d\n", 
                thisIndex,
                CkMyPe(),
                myNumParticles,
                myNumBuckets,
                totalNumTraversals,
                localTraversalState.pending,
                remoteTraversalState.pending);

  CkCallback cb(CkIndex_Main::quiescenceExit(),mainProxy);
  contribute(0,0,CkReduction::sum_int,cb);
}

void TreePiece::requestMoments(Key k, int replyTo){
  myDM->requestMoments(k,replyTo);
}

void TreePiece::requestParticles(RequestMsg *msg){
  // forward the request to the DM, since TPs don't
  // really own particles or nodes
  myDM->requestParticles(msg);
}

void TreePiece::requestNode(RequestMsg *msg){
  myDM->requestNode(msg);
}

int TreePiece::getIteration() {
  return iteration;
}

void TreePiece::pup(PUP::er &p){
  CBase_TreePiece::pup(p);
  p | iteration;
  if(p.isUnpacking()){
    init();
  }
}

#include "Traversal_defs.h"

