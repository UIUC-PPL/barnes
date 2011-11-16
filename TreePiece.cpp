#include "TreePiece.h"
#include "Messages.h"
#include "DataManager.h"
#include "Parameters.h"

#include "TaggedVector3D.h"
#include "Orb3dLB_notopo.h"

#include <fstream>

extern CProxy_Main mainProxy;
extern CProxy_DataManager dataManagerProxy;
extern Parameters globalParams;

TreePiece::TreePiece() {
  usesAtSync = CmiTrue;
  iteration = 0;
  findOrbLB();
  init();
}

void TreePiece::findOrbLB(){
  LBDatabase *lbdb = LBDatabaseObj();
  numLB = lbdb->getNLoadBalancers();
  BaseLB **lb = lbdb->getLoadBalancers();
  haveOrbLB = false;
  for(int i = 0; i < numLB; i++){
    if(string(lb[i]->lbName()) == "Orb3dLB_notopo"){
      orbLBProxy = lb[i]->getGroupID();
      haveOrbLB = true;
      break;
    }
  }
}

void TreePiece::init(){
  decompMsgsRecvd.length() = 0;
  totalNumTraversals = 2;
  myDM = dataManagerProxy.ckLocalBranch();
  root = NULL;
  myBuckets = NULL;
  myNumBuckets = 0;
  myNumParticles = 0;
  myRoot = NULL;
  localTraversalState.finishedIteration();
  remoteTraversalState.finishedIteration();
}

void TreePiece::receiveParticles(ParticleMsg *msg){
  int msgNumParticles = msg->numParticles;
  decompMsgsRecvd.push_back(msg);
  myNumParticles += msgNumParticles;
}

/*
void TreePiece::submitParticles(){
  if(thisIndex == 0) CkStartQD(CkCallback(CkIndex_TreePiece::submitParticles(),thisProxy));
  myDM->submitParticles(&decompMsgsRecvd,myNumParticles,this);
}
*/

void TreePiece::prepare(Node<ForceData> *_root, Node<ForceData> *_myRoot, Node<ForceData> **buckets, int _numBuckets){
  root = _root;
  myRoot = _myRoot;
  myBuckets = buckets;
  myNumBuckets = _numBuckets;
}

void TreePiece::startTraversal(int dummy){
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
    if(remoteTraversalState.complete()){
      doneRemoteRequests();
      traversalDone();
    }
  }
}

void TreePiece::doneRemoteRequests(){
  myDM->doneRemoteRequests();
  //dataManagerProxy[CkMyPe()].doneRemoteRequests();
}

void TreePiece::traversalDone(){ 
  totalNumTraversals--;
  if(totalNumTraversals == 0) finishIteration();
}

void TreePiece::finishIteration(){
#if 0
  Node<ForceData> **bucketptr;
  for(bucketptr = myBuckets; bucketptr != myBuckets+myNumBuckets; bucketptr++){
    Particle *p;
    ostringstream oss;
    for(p = (*bucketptr)->getParticles(); p != (*bucketptr)->getParticles()+(*bucketptr)->getNumParticles(); p++){
      oss << p->acceleration.x << " " << p->acceleration.y << " " << p->acceleration.z << " ; ";
    }
    CkPrintf("(%d,%d) bucket %d xfinal acc : %s\n", thisIndex, iteration, (*bucketptr)->getKey(), oss.str().c_str());
  }
#endif
  CmiUInt8 pn = localTraversalState.numInteractions[0]+remoteTraversalState.numInteractions[0];
  CmiUInt8 pp = localTraversalState.numInteractions[1]+remoteTraversalState.numInteractions[1];
  CmiUInt8 oc = localTraversalState.numInteractions[2]+remoteTraversalState.numInteractions[2];
  dataManagerProxy[CkMyPe()].traversalsDone(pn,pp,oc);

}

void TreePiece::cleanup(){
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

void TreePiece::requestParticles(std::pair<Key, int> request) {
  // forward the request to the DM, since TPs don't
  // really own particles or nodes
  myDM->requestParticles(request);
}

//void TreePiece::requestNode(RequestMsg *request){
void TreePiece::requestNode(std::pair<Key, int> &request) {
  myDM->requestNode(request);
}

int TreePiece::getIteration() {
  return iteration;
}

void TreePiece::pup(PUP::er &p){
  CBase_TreePiece::pup(p);
  p | iteration;
  if(p.isUnpacking()){
    init();
    findOrbLB();
  }
}

void TreePiece::startlb(){
  LDObjHandle handle = myRec->getLdHandle();
  if(numLB == 0){
    ResumeFromSync();
  }
  else if(haveOrbLB){
    Vector3D<float> centroid(0.0);
    if(myRoot != NULL) centroid = myRoot->data.moments.cm;

#if 0
    CkPrintf("tree piece %d pe %d contributing %f %f %f\n", 
                                        thisIndex,
                                        CkMyPe(),
                                        centroid.x,
                                        centroid.y,
                                        centroid.z
                                        );
#endif
    unsigned int ni = localTraversalState.numInteractions[0]
                     +localTraversalState.numInteractions[1]
                     +remoteTraversalState.numInteractions[0]
                     +remoteTraversalState.numInteractions[1];

    TaggedVector3D tv(centroid,handle,ni,myNumParticles,0,0);
    tv.tag = thisIndex;
    CkCallback lbcb(CkIndex_Orb3dLB_notopo::receiveCentroids(NULL),0,orbLBProxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else{
    // must call cleanup() before AtSync()
    cleanup();
    AtSync();
    return;
  }
  
  cleanup();
}

void TreePiece::doAtSync(){
  AtSync();
}

void TreePiece::ResumeFromSync() {
  //if(thisIndex == 0) CkPrintf("tree piece %d ResumeFromSync\n", thisIndex);

  CkCallback cb(CkIndex_DataManager::resumeFromLB(),dataManagerProxy);
  contribute(0,0,CkReduction::sum_int,cb);
}
#include "Traversal_defs.h"

