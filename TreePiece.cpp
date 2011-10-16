#include "TreePiece.h"
#include "Messages.h"
#include "DataManager.h"
#include "Parameters.h"

#include <fstream>

extern CProxy_Main mainProxy;
extern CProxy_DataManager dataManagerProxy;
extern Parameters globalParams;

TreePiece::TreePiece() : 
  myNumParticles(0),
  numDecompMsgsRecvd(0),
  smallestKey(~Key(0)),
  largestKey(Key(0)),
  localTraversalState(thisIndex,"local",this),
  remoteTraversalState(thisIndex,"remote",this),
  localStateID(0), 
  remoteStateID(1),
  localTraversalWorker(this),
  remoteTraversalWorker(this),
  totalNumTraversals(2),
  iteration(0),
  numTraversalsDone(0)
{
  myDM = dataManagerProxy.ckLocalBranch();
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

void TreePiece::prepare(Node<ForceData> *_root, Node<ForceData> **buckets, int bucketStart, int bucketEnd){
  root = _root;
  myBuckets = buckets+bucketStart;
  myNumBuckets = bucketEnd-bucketStart;
}

void TreePiece::startTraversal(){
  numTraversalsDone = 0;
  trav.setDataManager(myDM);

  if(myNumBuckets == 0){
    localGravityDone();
    remoteGravityDone();
    return;
  }

  remoteTraversalState.reset(myNumBuckets,myBuckets);
  remoteTraversalWorker.reset(&remoteTraversalState,*myBuckets);
  RescheduleMsg *msg = new (NUM_PRIORITY_BITS) RescheduleMsg;
  *(int *)CkPriorityPtr(msg) = REMOTE_GRAVITY_PRIORITY;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  thisProxy[thisIndex].doRemoteGravity(msg);

  localTraversalState.reset(myNumBuckets,myBuckets);
  localTraversalWorker.reset(&localTraversalState,*myBuckets);

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

  if(localTraversalState.decrPending(i)){
    localGravityDone();
    delete msg;
  }
  else{
    thisProxy[thisIndex].doLocalGravity(msg);
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

  if(remoteTraversalState.decrPending(i)){
    remoteGravityDone();
    delete msg;
  }
  else{
    thisProxy[thisIndex].doRemoteGravity(msg);
  }
}

void TreePiece::localGravityDone(){
  //CkPrintf("[%d] iter %d done local gravity state bucket %d/%d\n",thisIndex,getIteration(),localTraversalState.current,myNumBuckets);
  traversalDone();
}

void TreePiece::remoteGravityDone(){ 
  //CkPrintf("[%d] iter %d done local gravity state bucket %d/%d\n",thisIndex,getIteration(),remoteTraversalState.current,myNumBuckets);
  traversalDone();
}

void TreePiece::traversalDone(){ 
  numTraversalsDone++;

  if(numTraversalsDone == totalNumTraversals){
#ifdef STATISTICS
    CmiUInt8 pn = localTraversalState.numInteractions[0]+remoteTraversalState.numInteractions[0];
    CmiUInt8 pp = localTraversalState.numInteractions[1]+remoteTraversalState.numInteractions[1];
    CmiUInt8 oc = localTraversalState.numInteractions[2]+remoteTraversalState.numInteractions[2];
#ifdef CHECK_NUM_INTERACTIONS
    map<Key,CmiUInt8>::iterator it;
    for(it = localTraversalState.bucketNodeInteractions.begin(); 
        it != localTraversalState.bucketNodeInteractions.end();
        it++){
      myDM->addBucketNodeInteractions(it->first,it->second);
      //oss << "L " << it->first << " N inter " << it->second << endl;
    }
    for(it = localTraversalState.bucketPartInteractions.begin(); 
        it != localTraversalState.bucketPartInteractions.end();
        it++){
      myDM->addBucketPartInteractions(it->first,it->second);
      //oss << "L " << it->first << " P inter " << it->second << endl;
    }
    for(it = remoteTraversalState.bucketNodeInteractions.begin(); 
        it != remoteTraversalState.bucketNodeInteractions.end();
        it++){
      myDM->addBucketNodeInteractions(it->first,it->second);
      //oss << "R " << it->first << " N inter " << it->second << endl;
    }
    for(it = remoteTraversalState.bucketPartInteractions.begin(); 
        it != remoteTraversalState.bucketPartInteractions.end();
        it++){
      myDM->addBucketPartInteractions(it->first,it->second);
      //oss << "R " << it->first << " P inter " << it->second << endl;
    }
    CkPrintf("[%d] iter %d local %lu remote %lu\n", thisIndex, iteration, pn, pp);
#endif
    dataManagerProxy[CkMyPe()].traversalsDone(pn,pp,oc);
#else
    dataManagerProxy[CkMyPe()].traversalsDone();
#endif

    finishIteration();
  }
}

void TreePiece::finishIteration(){
  numTraversalsDone = 0;

  localTraversalState.finishedIteration();
  remoteTraversalState.finishedIteration();

  myNumParticles = 0;
  numDecompMsgsRecvd = 0;
  decompMsgsRecvd.length() = 0;

  iteration++;

  checkTraversals();
}

void TreePiece::checkTraversals(){
#ifdef DEBUG_TRAVERSALS
  map<Key,set<Key> >::iterator it;
  for(it = localTraversalState.bucketKeys.begin(); it != localTraversalState.bucketKeys.end(); it++){
    set<Key>::iterator iit;
    set<Key> &remlist = it->second;
    if(remlist.size() != 1){
      ostringstream oss;
      for(iit = remlist.begin(); iit != remlist.end(); iit++){
        oss << *iit << ",";
      }

      CkPrintf("[%d] local rem bucket %lu size %d: %s\n",
          thisIndex, it->first, remlist.size(), oss.str().c_str());
    }
    iit = remlist.find(Key(0));
    CkAssert(iit != remlist.end());
  }

  for(it = remoteTraversalState.bucketKeys.begin(); it != remoteTraversalState.bucketKeys.end(); it++){
    set<Key>::iterator iit;
    set<Key> &remlist = it->second;
    if(remlist.size() != 1){
      ostringstream oss;
      for(iit = remlist.begin(); iit != remlist.end(); iit++){
        oss << *iit << ",";
      }

      CkPrintf("[%d] remote rem bucket %lu size %d: %s\n",
          thisIndex, it->first, remlist.size(), oss.str().c_str());
    }
    iit = remlist.find(Key(0));
    CkAssert(iit != remlist.end());
  }
#endif
}

void TreePiece::quiescence(){
  CkPrintf("QUIESCENCE tree piece %d proc %d submitted %d numBuckets %d trav_done %d outstanding local %d remote %d\n", 
                thisIndex,
                CkMyPe(),
                myNumParticles,
                myNumBuckets,
                numTraversalsDone,
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

#include "Traversal_defs.h"

