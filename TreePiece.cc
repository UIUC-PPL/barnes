/*
 * CharmBH: TreePiece.cc
 * TreePieces perform parallel traversals on the 
 * distributed tree data structure. They do not
 * make direct requests for data, but instead 
 * rely on Traversal objects to do so.
 *
 * A TreePiece traverses the global tree once for 
 * each of its buckets. The traversal begins at the
 * root, and for each node encountered, the traversal
 * checks whether the node is far enough away from the
 * bucket that the interaction of its particles and 
 * those of the bucket can be approximated through a 
 * bucket node interaction. If so, the computation
 * is performed. If not, the traversal attempts to 
 * consider each of the children of the node in turn.
 * However, if the node turns out to be a leaf, 
 * pair-wise interactions are performed between the
 * particles of the node and those of the bucket.
 */
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
  myRoot = NULL;
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
  root = NULL;
  myBuckets = NULL;
  myNumBuckets = 0;
}

// Particles are decomposed onto TreePieces, and the following
// method is invoked to communicate the particles to them
void TreePiece::receiveParticles(ParticleMsg *msg){
  int msgNumParticles = msg->numParticles;
  decompMsgsRecvd.push_back(msg);
  myNumParticles += msgNumParticles;
  numDecompMsgsRecvd++;
  if(smallestKey > msg->part[0].key) smallestKey = msg->part[0].key;
  if(largestKey < msg->part[msgNumParticles-1].key) largestKey = msg->part[msgNumParticles-1].key;
  
  // If I have received all the messages I need, I can submit 
  // my particles to the DataManager on this PE
  if(numDecompMsgsRecvd == CkNumPes()){
    submitParticles();
    numDecompMsgsRecvd = 0;
  }
}

// If a PE has no particles for me, it invokes this method
void TreePiece::receiveParticles(){
  numDecompMsgsRecvd++;
  if(numDecompMsgsRecvd == CkNumPes()){
    submitParticles();
    numDecompMsgsRecvd = 0;
  }
}

// Submit particles to data manager
void TreePiece::submitParticles(){
  myDM->submitParticles(&decompMsgsRecvd,myNumParticles,this,smallestKey,largestKey);
}

// Invoked on me by the DataManager on my PE
void TreePiece::prepare(Node<ForceData> *_root, Node<ForceData> *_myRoot, Node<ForceData> **buckets, int bucketStart, int bucketEnd){
  root = _root;
  myRoot = _myRoot;
  myBuckets = buckets+bucketStart;
  myNumBuckets = bucketEnd-bucketStart;
}

// Invoked on me by my DataManager 
void TreePiece::startTraversal(){
  // Initialize the Traversal object
  trav.setDataManager(myDM);
  
  // If I received no particles, I can finish iteration
  if(myNumBuckets == 0){
    finishIteration();
    return;
  }

  // Set up the remote traversal: this traversal only processes
  // Remote nodes
  remoteTraversalState.reset(this,myNumBuckets,myBuckets);
  remoteTraversalWorker.reset(this,&remoteTraversalState,*myBuckets);
  RescheduleMsg *msg = new (NUM_PRIORITY_BITS) RescheduleMsg;
  *(int *)CkPriorityPtr(msg) = REMOTE_GRAVITY_PRIORITY;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  thisProxy[thisIndex].doRemoteGravity(msg);

  // Set up local traversal: this traversal only processes 
  // Local nodes
  localTraversalState.reset(this,myNumBuckets,myBuckets);
  localTraversalWorker.reset(this,&localTraversalState,*myBuckets);
  msg = new (NUM_PRIORITY_BITS) RescheduleMsg;
  *(int *)CkPriorityPtr(msg) = LOCAL_GRAVITY_PRIORITY;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  thisProxy[thisIndex].doLocalGravity(msg);
}

// Self-scheduling local gravity computation. Relinquishes processor
// every yieldPeriod buckets.
void TreePiece::doLocalGravity(RescheduleMsg *msg){
  int i;
  // Iterate through buckets
  for(i = 0; i < globalParams.yieldPeriod && 
                 localTraversalState.current < myNumBuckets; 
                 i++){
    // Initiate a top-down traversal of the global tree (rooted at 'root')
    // with the LocalTraversalWorker. The book-keeping state for this traversal
    // is kept in localTraversalState.
    trav.topDownTraversal(root,&localTraversalWorker,&localTraversalState);
    // Finished traversal for this bucket, move to next one.
    localTraversalState.current++;
    localTraversalState.currentBucketPtr++;
    // Tell the worker to move to the next bucket
    localTraversalWorker.setContext(*localTraversalState.currentBucketPtr);
  }

  // We are one bucket closer to completing the traversal
  localTraversalState.decrPending(i);
  // Are we done with all the buckets? 
  if(localTraversalState.current < myNumBuckets){
    CkAssert(!localTraversalState.complete());
    // If not, schedule a few more buckets for local traversal
    thisProxy[thisIndex].doLocalGravity(msg);
  }
  else{
    delete msg;
    // Traversal complete
    if(localTraversalState.complete()) traversalDone();
  }
}

// Ditto for Remote traversal
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
    // It is unlikely that the Remote traversal will have
    // completed here, because it most probably triggered
    // some remote communication requests.
    if(remoteTraversalState.complete()) traversalDone();
  }
}

// If we are done with all (both) of our traversals, we are
// done with the iteration 
void TreePiece::traversalDone(){ 
  totalNumTraversals--;
  if(totalNumTraversals == 0) finishIteration();
}

// Tell the DataManager on this PE that I am done with my traversal
void TreePiece::finishIteration(){
  CmiUInt8 pn = localTraversalState.numInteractions[0]+remoteTraversalState.numInteractions[0];
  CmiUInt8 pp = localTraversalState.numInteractions[1]+remoteTraversalState.numInteractions[1];
  CmiUInt8 oc = localTraversalState.numInteractions[2]+remoteTraversalState.numInteractions[2];
  dataManagerProxy[CkMyPe()].traversalsDone(pn,pp,oc);

  // Reset local and remote traversal book-keeping states.
  localTraversalState.finishedIteration();
  remoteTraversalState.finishedIteration();
  
  init();

  iteration++;
}

// These functions just pass on requests for data assigned
// to me to my DataManager, since the DM keeps track of all
// data on the PE.
void TreePiece::requestMoments(Key k, int replyTo){
  myDM->requestMoments(k,replyTo);
}

void TreePiece::requestParticles(std::pair<Key, int> request) {
  // forward the request to the DM, since TPs don't
  // really own particles or nodes
  myDM->requestParticles(request);
}

void TreePiece::requestNode(std::pair<Key, int> request) {
  myDM->requestNode(request);
}

// A TreePiece doesn't have much state to save, especially
// since load balancing is only done at iteration boundaries.
void TreePiece::pup(PUP::er &p){
  CBase_TreePiece::pup(p);
  p | iteration;
  if(p.isUnpacking()){
    init();
  }
}

// Commence load balancing
void TreePiece::startlb(){
  AtSync();
}

// Resume after load balancing: transfer control to data manager
void TreePiece::ResumeFromSync() {
  CkCallback cb(CkIndex_DataManager::resumeFromLB(),dataManagerProxy);
  contribute(0,0,CkReduction::sum_int,cb);
}
#include "Traversal_defs.h"

