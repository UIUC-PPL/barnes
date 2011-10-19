#include "DataManager.h"
#include "Reduction.h"
#include "defines.h"
#include "Messages.h"
#include "Parameters.h"

#include "Worker.h"
#include "TreePiece.h"

#include "Request.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
extern CProxy_TreePiece treePieceProxy;
extern CProxy_Main mainProxy;
extern Parameters globalParams;

void copyMomentsToNode(Node<ForceData> *node, const MomentsExchangeStruct &mes){
  CkAssert(node->getKey() == mes.key);

  node->data.moments = mes.moments;
  node->data.box = mes.box;
  NodeType type = mes.type;
  node->setType(Node<ForceData>::makeRemote(type));

}

DataManager::DataManager() : 
  iteration(0),
  prevIterationStart(0.0),
  savedEnergy(0.0),
  root(NULL),
  keyRanges(NULL), 
  rangeMsg(NULL)
{
  init();
}

void DataManager::loadParticles(CkCallback &cb){
  numRankBits = LOG_BRANCH_FACTOR;

  const char *fname = globalParams.filename.c_str();
  int npart = globalParams.numParticles;

  std::ifstream partFile;
  partFile.open(fname, ios::in | ios::binary);
  CkAssert(partFile.is_open());

  int offset = 0; 
  int myid = CkMyPe();
  int npes = CkNumPes();

  int avgParticlesPerPE = npart/npes;
  int rem = npart-npes*avgParticlesPerPE;
  if(myid < rem){
    avgParticlesPerPE++;
    offset = myid*avgParticlesPerPE;
  }
  else{
    offset = myid*avgParticlesPerPE+rem;
  }
  myNumParticles = avgParticlesPerPE;
  offset *= SIZE_PER_PARTICLE;
  offset += PREAMBLE_SIZE;

  myParticles.reserve(myNumParticles);
  myParticles.length() = myNumParticles;

  partFile.clear();
  partFile.seekg(offset,ios::beg);
  if(partFile.fail()){
    std::ostringstream oss;
    oss << "couldn't seek to position " << offset << " on PE " << CkMyPe() << " position " << partFile.tellg() << endl;
    CkAbort(oss.str().c_str());
  }
  unsigned int numParticlesDone = 0;

  BoundingBox myBox;

  Real tmp[REALS_PER_PARTICLE];
  while(numParticlesDone < myNumParticles && !partFile.eof()){
    partFile.read((char *)tmp, SIZE_PER_PARTICLE);
    myParticles[numParticlesDone].position.x = tmp[0];
    myParticles[numParticlesDone].position.y = tmp[1];
    myParticles[numParticlesDone].position.z = tmp[2];
    myParticles[numParticlesDone].velocity.x = tmp[3];
    myParticles[numParticlesDone].velocity.y = tmp[4];
    myParticles[numParticlesDone].velocity.z = tmp[5];
    myParticles[numParticlesDone].mass = tmp[6];

    myParticles[numParticlesDone].acceleration.x = 0.0;
    myParticles[numParticlesDone].acceleration.y = 0.0;
    myParticles[numParticlesDone].acceleration.z = 0.0;
    myParticles[numParticlesDone].potential = 0.0;
    myBox.grow(myParticles[numParticlesDone].position);

    numParticlesDone++;
  }

  CkAssert(numParticlesDone == myNumParticles);
  myBox.numParticles = myNumParticles;

  partFile.close();

  contribute(sizeof(BoundingBox),&myBox,BoundingBoxGrowReductionType,cb);
}

void DataManager::hashParticleCoordinates(OrientedBox<Real> &universe){
  Key prepend;
  prepend = 1L;
  prepend <<= (TREE_KEY_BITS-1);

  Real xsz = universe.greater_corner.x-universe.lesser_corner.x;
  Real ysz = universe.greater_corner.y-universe.lesser_corner.y;
  Real zsz = universe.greater_corner.z-universe.lesser_corner.z;

  for(unsigned int i = 0; i < myNumParticles; i++){
    Particle *p = &(myParticles[i]);
    Key xint = ((Key) (((p->position.x-universe.lesser_corner.x)*(BOXES_PER_DIM*1.0))/xsz)); 
    Key yint = ((Key) (((p->position.y-universe.lesser_corner.y)*(BOXES_PER_DIM*1.0))/ysz)); 
    Key zint = ((Key) (((p->position.z-universe.lesser_corner.z)*(BOXES_PER_DIM*1.0))/zsz)); 

    Key mask = Key(0x1);
    Key k = Key(0x0);
    int shiftBy = 0;
    for(int j = 0; j < BITS_PER_DIM; j++){
      k |= ((zint & mask) <<  shiftBy);
      k |= ((yint & mask) << (shiftBy+1));
      k |= ((xint & mask) << (shiftBy+2));
      mask <<= 1;
      // minus 1 because mask itself has shifted
      // left by one position
      shiftBy += (NDIMS-1);
    }
    k |= prepend; 
    myParticles[i].key = k;
  }
}

void DataManager::decompose(BoundingBox &universe){

  hashParticleCoordinates(universe.box);
  myParticles.quickSort();

  if(CkMyPe()==0){
    float memMB = (1.0*CmiMemoryUsage())/(1<<20);
    ostringstream oss; 
    CkPrintf("(%d) prev time %g s\n", CkMyPe(), CmiWallTimer()-prevIterationStart);
    CkPrintf("(%d) start iteration %d\n", CkMyPe(), iteration);
    CkPrintf("(%d) mem %.2f MB\n", CkMyPe(), memMB);
    CkPrintf("(%d) univ %f %f %f %f %f %f energy %f\n", 
              CkMyPe(),
              universe.box.lesser_corner.x,
              universe.box.lesser_corner.y,
              universe.box.lesser_corner.z,
              universe.box.greater_corner.x,
              universe.box.greater_corner.y,
              universe.box.greater_corner.z,
              universe.energy);
 

    prevIterationStart = CkWallTimer();
  }

  numTreePieces = 1;
  initHistogramParticles();
  sendHistogram();
}

void DataManager::initHistogramParticles(){
  int rootDepth = 0;
  
  sortingRoot = new Node<NodeDescriptor>(Key(1),
                         rootDepth,
                         myParticles.getVec(),
                         myNumParticles);
  activeBins.addNewNode(sortingRoot);

  // don't access myParticles through ckvec after this
  // anyway. these must be reset before this DM starts
  // to receive submitted particles from TPs placed on it
  myNumParticles = 0;
  myParticles.length() = 0;
}

void DataManager::sendHistogram(){
  CkCallback cb(CkIndex_DataManager::receiveHistogram(NULL),0,this->thisgroup);
  contribute(sizeof(NodeDescriptor)*activeBins.getNumCounts(),activeBins.getCounts(),NodeDescriptorReductionType,cb);
  activeBins.reset();
}

// executed on PE 0
void DataManager::receiveHistogram(CkReductionMsg *msg){
  int numRecvdBins = msg->getSize()/sizeof(NodeDescriptor);
  NodeDescriptor *descriptors = (NodeDescriptor *)msg->getData();
  CkVec<int> binsToRefine;

  binsToRefine.reserve(numRecvdBins);
  binsToRefine.length() = 0;

  int particlesHistogrammed = 0;

  CkVec<pair<Node<NodeDescriptor>*,bool> > *active = activeBins.getActive();
  CkAssert(numRecvdBins == active->length());

  for(int i = 0; i < numRecvdBins; i++){
    if(descriptors[i].numParticles > (Real)(DECOMP_TOLERANCE*globalParams.ppc)){
      // need to refine this bin (partition)
      binsToRefine.push_back(i);
      numTreePieces += (BRANCH_FACTOR-1);
      if(numTreePieces > globalParams.numTreePieces){
        CkPrintf("have %d treepieces need %d\n",globalParams.numTreePieces,numTreePieces);
        CkAbort("Need more tree pieces!\n");
      }
    }
    else{
      Node<NodeDescriptor> *nd = (*active)[i].first;
      nd->data = descriptors[i];
    }

    particlesHistogrammed += descriptors[i].numParticles;
  }

  if(binsToRefine.size()) {
    thisProxy.receiveSplitters(binsToRefine);
    decompIterations++;
  }
  else{
    // create tree pieces and send proxy
    CkPrintf("[0] decomp done after %d iterations used treepieces %d\n", decompIterations, numTreePieces);
    decompIterations = 0;
    
    //keyRanges.reserve(numTreePieces*2);
    keyRanges = new Key[numTreePieces*2];
    // so that by the time tree pieces start submitting
    // particles (which can only happen after the flushParticles()
    // below, we have the right count of local tree pieces
    senseTreePieces();
    flushParticles();
    // PE 0 sets ranges in sendParticlesToTreePiece
    haveRanges = true;

    int numKeys = numTreePieces*2;

    RangeMsg *rmsg = new (numKeys) RangeMsg;
    rmsg->numTreePieces = numTreePieces;
    memcpy(rmsg->keys,keyRanges,sizeof(Key)*numKeys);
    thisProxy.sendParticles(rmsg);
  }

  delete msg;
}

void DataManager::flushParticles(){
  ParticleFlushWorker pfw(this);
  scaffoldTrav.preorderTraversal(sortingRoot,&pfw);

  int numUsefulTreePieces = pfw.getNumLeaves(); 
  for(int i = numUsefulTreePieces; i < globalParams.numTreePieces; i++){
    treePieceProxy[i].receiveParticles();
  }

  // done with sorting tree; delete
  FreeTreeWorker<NodeDescriptor> freeWorker;
  scaffoldTrav.postorderTraversal(sortingRoot,&freeWorker);

  delete sortingRoot;
  sortingRoot = NULL;

}

void DataManager::receiveSplitters(CkVec<int> splitBins) {

  // process bins to refine
  activeBins.processRefine(splitBins.getVec(), splitBins.size());

  // We traverse the final tree to flush particles to 
  // appropriate tree pieces

  sendHistogram();
}

void DataManager::sendParticlesToTreePiece(Node<NodeDescriptor> *nd, int tp) {
  CkAssert(nd->getNumChildren() == 0);
  int np = nd->getNumParticles();

  if(np > 0){
    ParticleMsg *msg = new (np,0) ParticleMsg;
    memcpy(msg->part, nd->getParticles(), sizeof(Particle)*np);
    msg->numParticles = np;
    treePieceProxy[tp].receiveParticles(msg);
  }
  else{
    treePieceProxy[tp].receiveParticles();
  }

  // only PE 0 has the correct ranges
  if(CkMyPe() == 0){
    if(nd->data.numParticles > 0){
      CkAssert(nd->data.smallestKey <= nd->data.largestKey);
    } else {
      CkAssert(nd->data.smallestKey == nd->data.largestKey);
    }

    keyRanges[(tp<<1)] = nd->data.smallestKey;
    keyRanges[(tp<<1)+1] = nd->data.largestKey;

  }
}

void DataManager::sendParticles(RangeMsg *msg){
  if(CkMyPe() != 0){
    numTreePieces = msg->numTreePieces;
    keyRanges = msg->keys;
    haveRanges = true;
    rangeMsg = msg;
    flushParticles();

    senseTreePieces();

    if(submittedParticles.length() == numLocalTreePieces){
      processSubmittedParticles();
    }
  }
  else{
    CkAssert(numTreePieces == msg->numTreePieces);
    CkAssert(haveRanges);
    delete msg;
  }
}

void DataManager::senseTreePieces(){
  localTreePieces.reset();
  CkLocMgr *mgr = treePieceProxy.ckLocMgr();
  mgr->iterate(localTreePieces);
  numLocalTreePieces = localTreePieces.count;
}

void DataManager::submitParticles(CkVec<ParticleMsg*> *vec, int numParticles, TreePiece * tp, Key smallestKey, Key largestKey){ 
  submittedParticles.push_back(TreePieceDescriptor(vec,numParticles,tp,tp->getIndex(),smallestKey,largestKey));
  myNumParticles += numParticles;
  if(submittedParticles.length() == numLocalTreePieces && haveRanges){
    processSubmittedParticles();
  }
}

void DataManager::processSubmittedParticles(){
  int offset = 0;

  submittedParticles.quickSort();
  
  myParticles.resize(myNumParticles);

  for(int i = 0; i < submittedParticles.length(); i++){
    TreePieceDescriptor &descr = submittedParticles[i];
    CkVec<ParticleMsg*> *vec = descr.vec;
    for(int j = 0; j < vec->length(); j++){
      ParticleMsg *msg = (*vec)[j];
      memcpy(myParticles.getVec()+offset,msg->part,sizeof(Particle)*msg->numParticles);
      offset += msg->numParticles;
      delete msg;
    }
  }

  myParticles.quickSort();

  buildTree();
  // add dummy tree piece whose index is larger than
  // that of all others. this is required to mark the
  // boundary of nodes/particles owned by this PE.
  submittedParticles.push_back(TreePieceDescriptor(globalParams.numTreePieces));

  // makeMoments also sends out requests for moments 
  // of remote nodes
  makeMoments();

#if 0
  ostringstream oss;
  oss << "tree." << CkMyPe() << ".dot";

  ofstream ofs(oss.str().c_str());
  ofs << "digraph PE" << CkMyPe() << "{" << endl;
  if(root != NULL) printTree(root,ofs);
  ofs << "}" << endl;
  ofs.close();
#endif

  doneTreeBuild = true;

  // are all particles local to this PE? 
  if(root != NULL && root->getType() == Internal){
    passMomentsUpward(root);
  }

  flushMomentRequests();
}

void DataManager::buildTree(){

  int rootDepth = 0;
  root = new Node<ForceData>(Key(1),rootDepth,myParticles.getVec(),myNumParticles);
  root->setOwners(0,numTreePieces-1);
  nodeTable[Key(1)] = root;
  if(myNumParticles == 0){
    return;
  }

  OwnershipActiveBinInfo<ForceData> abi(keyRanges);
  abi.addNewNode(root);
  int numFatNodes = 1;

  int limit = ((Real)globalParams.ppb*BUCKET_TOLERANCE);

  CkVec<int> refines;

  while(numFatNodes > 0){
    abi.reset();
    refines.length() = 0;

    CkVec<std::pair<Node<ForceData>*,bool> > *active = abi.getActive();

    for(int i = 0; i < active->length(); i++){
      Node<ForceData> *node = (*active)[i].first;
      if(node->getNumParticles() > 0 && ((node->getOwnerEnd() > node->getOwnerStart()) || 
         (node->getNumParticles() > limit))){
        refines.push_back(i);
      }
    }

    abi.processRefine(refines.getVec(), refines.length());
    numFatNodes = abi.getNumCounts();
  }

}

void DataManager::makeMoments(){
  if(root == NULL) return;

  MomentsWorker mw(submittedParticles, nodeTable, localTPRoots, myBuckets);
  fillTrav.postorderTraversal(root,&mw);
}

Node<ForceData> *DataManager::lookupNode(Key k){
  map<Key,Node<ForceData>*>::iterator it;
  it = nodeTable.find(k);
  if(it == nodeTable.end()) return NULL;
  else return it->second;
}

void DataManager::requestMoments(Key k, int replyTo){
  pendingMoments[k].push_back(replyTo);
  TB_DEBUG("(%d) received requestMoments %lu from pe %d doneTreeBuild %d\n", 
          CkMyPe(), k, replyTo, doneTreeBuild);

  if(doneTreeBuild){    
    Node<ForceData> *node = lookupNode(k);
    if(node == NULL){
      CkPrintf("(%d) recvd request from %d for moments %lu\n", CkMyPe(), replyTo, k);
      CkAbort("bad request\n");
    }
    bool ready = node->allChildrenMomentsReady();
    TB_DEBUG("(%d) node %lu ready %d\n", CkMyPe(), k, ready);

    if(ready){
      map<Key,CkVec<int> >::iterator it = pendingMoments.find(k);
      CkVec<int> &requestors = it->second;
      respondToMomentsRequest(node,requestors);
      pendingMoments.erase(it);
    }
  }
}

void DataManager::flushMomentRequests(){
  CkAssert(doneTreeBuild);
  map<Key,CkVec<int> >::iterator it;
  for(it = pendingMoments.begin(); it != pendingMoments.end();){
    Key k = it->first;
    Node<ForceData> *node = lookupNode(k);
    CkAssert(node != NULL);
    if(node->allChildrenMomentsReady()){
      CkVec<int> &requestors = it->second;
      respondToMomentsRequest(node,requestors);
      map<Key,CkVec<int> >::iterator kill = it;
      ++it;
      pendingMoments.erase(kill);
    }
    else{
      ++it;
    }
  }
}

void DataManager::respondToMomentsRequest(Node<ForceData> *node, CkVec<int> &replyTo){
  MomentsExchangeStruct moments = *node;
  for(int i = 0; i < replyTo.length(); i++){
    TB_DEBUG("(%d) responding to %d with node %lu\n", CkMyPe(), replyTo[i], node->getKey());
    thisProxy[replyTo[i]].receiveMoments(moments);
  }
  replyTo.length() = 0;
}

void DataManager::receiveMoments(MomentsExchangeStruct moments){
  Node<ForceData> *node = lookupNode(moments.key);
  CkAssert(node != NULL);

  // update moments of leaf and pass these on 
  // to parent recursively; if there are requests
  // for these nodes, respond to them
  updateLeafMoments(node,moments);
}

void DataManager::updateLeafMoments(Node<ForceData> *node, MomentsExchangeStruct &data){
  copyMomentsToNode(node,data);
  TB_DEBUG("(%d) updateLeafMoments %lu\n", CkMyPe(), node->getKey());
  passMomentsUpward(node);
}

void DataManager::passMomentsUpward(Node<ForceData> *node){
  TB_DEBUG("(%d) passUp %lu\n", CkMyPe(), node->getKey());
  map<Key,CkVec<int> >::iterator it = pendingMoments.find(node->getKey());
  if(it != pendingMoments.end()){
    CkVec<int> &requestors = it->second;
    respondToMomentsRequest(node,requestors);
    pendingMoments.erase(it);
  }

  Node<ForceData> *parent = node->getParent();
  if(parent == NULL){
    CkAssert(node->getKey() == Key(1));
    treeReady();
  }else{
    parent->childMomentsReady();
    TB_DEBUG("[%d] parent %lu children ready %d\n", CkMyPe(), parent->getKey(), parent->getNumChildrenMomentsReady());
    if(parent->allChildrenMomentsReady()){
      parent->getMomentsFromChildren();
      parent->getOwnershipFromChildren();
      passMomentsUpward(parent);
    }
  }
}

// doneTreeBuild: built local tree and sent out requests for remote 

void DataManager::treeReady(){
  treeMomentsReady = true;
  flushBufferedRemoteDataRequests();
  startTraversal();
}

void DataManager::flushBufferedRemoteDataRequests(){
  CkAssert(treeMomentsReady);
  for(int i = 0; i < bufferedNodeRequests.length(); i++){
    requestNode(bufferedNodeRequests[i]);
  }
  for(int i = 0; i < bufferedParticleRequests.length(); i++){
    requestParticles(bufferedParticleRequests[i]);
  }
  bufferedNodeRequests.length() = 0;
  bufferedParticleRequests.length() = 0;
}

bool CompareNodePtrToKey(void *a, Key k){
  Node<ForceData> *node = *((Node<ForceData>**)a);
  return (Node<ForceData>::getParticleLevelKey(node) >= k);
}

void DataManager::startTraversal(){
  Node<ForceData> **bucketPtrs = myBuckets.getVec();
  submittedParticles[0].bucketStartIdx = 0;
  int start = 0;
  int end = myBuckets.length();

  if(end > 0){
    for(int i = 0; i < numLocalTreePieces-1; i++){
      TreePieceDescriptor &descr = submittedParticles[i];
      int bucketIdx = binary_search_ge<Node<ForceData>*>(descr.largestKey,bucketPtrs,start,end,CompareNodePtrToKey);
      descr.bucketEndIdx = bucketIdx;
      int tpIndex = descr.owner->getIndex();
      descr.owner->prepare(root,localTPRoots[tpIndex],myBuckets.getVec(),descr.bucketStartIdx,descr.bucketEndIdx);
      treePieceProxy[tpIndex].startTraversal();
      submittedParticles[i+1].bucketStartIdx = bucketIdx;
      start = bucketIdx;
    }
    TreePieceDescriptor &descr = submittedParticles[numLocalTreePieces-1];
    descr.bucketEndIdx = myBuckets.length();
    int tpIndex = descr.owner->getIndex();
    descr.owner->prepare(root,localTPRoots[tpIndex],myBuckets.getVec(),descr.bucketStartIdx,descr.bucketEndIdx);
    treePieceProxy[tpIndex].startTraversal();
  }
  else if(numLocalTreePieces > 0){
    for(int i = 0; i < numLocalTreePieces; i++){
      TreePieceDescriptor &descr = submittedParticles[i];
      descr.owner->prepare(root,NULL,myBuckets.getVec(),0,0);
      int tpIndex = descr.owner->getIndex();
      treePieceProxy[tpIndex].startTraversal();
    }
  }
  else{
    finishIteration();
  }
}

void DataManager::requestParticles(Node<ForceData> *leaf, CutoffWorker<ForceData> *worker, State *state, Traversal<ForceData> *traversal){
  Key key = leaf->getKey();
  Request &request = particleRequestTable[key];
#ifdef TRACE_REMOTE_DATA_REQUESTS
  traceUserEvent(REMOTE_PARTICLE_REQUEST);
#endif
  if(!request.sent){
    partReqs.incrRequests();

    if(leaf->isCached()) request.parentCached = true;
    else request.parentCached = false;

    CkAssert(leaf->getOwnerStart() == leaf->getOwnerEnd());
    int owner = leaf->getOwnerStart();
    treePieceProxy[owner].requestParticles( std::make_pair(key, CkMyPe()) );
    request.sent = true;
    
    request.parent = leaf;
    RRDEBUG("(%d) REQUEST particles %lu from tp %d\n", CkMyPe(), key, owner);
  }
  request.requestors.push_back(Requestor(worker,state,traversal,worker->getContext()));
  partReqs.incrDeliveries();
}

void DataManager::requestParticles(std::pair<Key, int> request) {
  if(!treeMomentsReady){
    bufferedParticleRequests.push_back(request);
    return;
  }

  RRDEBUG("(%d) REPLY particles key %lu to %d\n", CkMyPe(), request.first, request.second);

  map<Key,Node<ForceData>*>::iterator it = nodeTable.find(request.first);
  CkAssert(it != nodeTable.end());
  Node<ForceData> *bucket = it->second;
  CkAssert(bucket->getType() == Bucket);

  Particle *data = bucket->getParticles();
  int np = bucket->getNumParticles();

  ParticleReplyMsg *pmsg = new (np,NUM_PRIORITY_BITS) ParticleReplyMsg;
  *(int *)CkPriorityPtr(pmsg) = RECV_PARTICLES_PRIORITY;
  CkSetQueueing(pmsg,CK_QUEUEING_IFIFO);

  pmsg->key = request.first;
  pmsg->np = np;
  for(int i = 0; i < np; i++){
    pmsg->data[i] = data[i];
  }

  thisProxy[request.second].recvParticles(pmsg);
}

void DataManager::requestNode(Node<ForceData> *leaf, CutoffWorker<ForceData> *worker, State *state, Traversal<ForceData> *traversal){
  Key key = leaf->getKey();
  Request &request = nodeRequestTable[key];
#ifdef TRACE_REMOTE_DATA_REQUESTS
  traceUserEvent(REMOTE_NODE_REQUEST);
#endif
  if(!request.sent){
    nodeReqs.incrRequests();

    if(leaf->isCached()) request.parentCached = true;
    else request.parentCached = false;

    int numOwners = leaf->getOwnerEnd()-leaf->getOwnerStart()+1;
    int requestOwner = leaf->getOwnerStart()+(rand()%numOwners);
    RRDEBUG("(%d) REQUEST node %lu from tp %d\n", CkMyPe(), key, requestOwner);
    treePieceProxy[requestOwner].requestNode( std::make_pair(key, CkMyPe()) );
    request.sent = true;
    request.parent = leaf;
  }
  request.requestors.push_back(Requestor(worker,state,traversal,worker->getContext()));
  nodeReqs.incrDeliveries();
}

void DataManager::requestNode(std::pair<Key, int> request) {
  if(!treeMomentsReady){
    bufferedNodeRequests.push_back(request);
    return;
  }

  RRDEBUG("(%d) REPLY node %lu to %d\n", CkMyPe(), request.first, request.second);


  map<Key,Node<ForceData>*>::iterator it = nodeTable.find(request.first);
  CkAssert(it != nodeTable.end());
  Node<ForceData> *node = it->second;

  CkAssert(node->getNumChildren() > 0);
#if 0
  if(node->getNumChildren() == 0){
    CkPrintf("[%d] children of leaf node %lu type %s requested!\n", CkMyPe(), node->getKey(), NodeTypeString[node->getType()].c_str());
    CkAbort("Leaf children request\n");
  }
#endif

  TreeSizeWorker tsz(node->getDepth()+globalParams.cacheLineSize);
  fillTrav.topDownTraversal_local(node,&tsz);

  int nn = tsz.getNumNodes();

  NodeReplyMsg *nmsg = new (nn,NUM_PRIORITY_BITS) NodeReplyMsg;
  *(int *)CkPriorityPtr(nmsg) = RECV_NODE_PRIORITY;
  CkSetQueueing(nmsg,CK_QUEUEING_IFIFO);

  nmsg->key = request.first;
  nmsg->nn = nn;

  Node<ForceData> *emptyBuf = nmsg->data;
  node->serialize(NULL,emptyBuf,globalParams.cacheLineSize);
  CkAssert(emptyBuf == nmsg->data+nn);

  thisProxy[request.second].recvNode(nmsg);
}

void DataManager::recvParticles(ParticleReplyMsg *msg){
  map<Key,Request>::iterator it = particleRequestTable.find(msg->key);
  CkAssert(it != particleRequestTable.end());

  RRDEBUG("(%d) RECVD particle REPLY for key %lu\n", CkMyPe(), msg->key);

  Request &req = it->second;
  CkAssert(req.requestors.length() > 0);
  CkAssert(req.sent);
  CkAssert(req.msg == NULL);

  req.msg = msg;
  req.data = msg->data;
  
  // attach particles to bucket in tree 
  Node<ForceData> *leaf = req.parent;
  CkAssert(leaf != NULL);
  CkAssert(leaf->getType() == RemoteBucket);
  leaf->setParticles((Particle *)msg->data,msg->np);

  partReqs.decrRequests();
  partReqs.decrDeliveries(req.requestors.length());
  req.deliverParticles(msg->np);
}

void DataManager::recvNode(NodeReplyMsg *msg){
  map<Key,Request>::iterator it = nodeRequestTable.find(msg->key);
  CkAssert(it != nodeRequestTable.end());

  RRDEBUG("(%d) RECVD node REPLY for key %lu\n", CkMyPe(), msg->key);

  Request &req = it->second;
  CkAssert(req.requestors.length() > 0);
  CkAssert(req.sent);
  CkAssert(req.msg == NULL);

  req.msg = msg;
  req.data = msg->data;
  
  // attach recvd subtree to appropriate point in local tree
  Node<ForceData> *node = req.parent;
  CkAssert(node != NULL);
  
  node->deserialize(msg->data, msg->nn);

  nodeReqs.decrRequests();
  nodeReqs.decrDeliveries(req.requestors.length());
  RRDEBUG("(%d) DELIVERING key %lu\n", CkMyPe(), msg->key);

  req.deliverNode();
  RRDEBUG("(%d) DELIVERED key %lu\n", CkMyPe(), msg->key);
}

void DataManager::traversalsDone(CmiUInt8 pnInter, CmiUInt8 ppInter, CmiUInt8 openCrit)
{
  numTreePiecesDoneTraversals++;
  numInteractions[0] += pnInter;
  numInteractions[1] += ppInter;
  numInteractions[2] += openCrit;
  if(numTreePiecesDoneTraversals == numLocalTreePieces){
    finishIteration();
  }
}

void DataManager::finishIteration(){
  // can't advance particles here, because other PEs 
  // might not have finished their traversals yet, 
  // and therefore might need my particles

  CkAssert(nodeReqs.test());
  CkAssert(partReqs.test());

  DtReductionStruct dtred;
  findMinVByA(dtred);

  dtred.pnInteractions = numInteractions[0];
  dtred.ppInteractions = numInteractions[1];
  dtred.openCrit = numInteractions[2];

  CkCallback cb(CkIndex_DataManager::advance(NULL),thisProxy);
  contribute(sizeof(DtReductionStruct),&dtred,DtReductionType,cb);

}

void DataManager::advance(CkReductionMsg *msg){

  DtReductionStruct *dtred = (DtReductionStruct *)(msg->getData());
  if(dtred->haveNaN){
    CkPrintf("(%d) iteration %d NaN accel detected! Exit...\n", CkMyPe(), iteration);
    markNaNBuckets();
    CkCallback exitCb(CkCallback::ckExit);
    contribute(0,0,CkReduction::sum_int,exitCb);
    return;
  }

  myBox.reset();
  kickDriftKick(myBox.box,myBox.energy);

  Real pad = 0.001;
  myBox.expand(pad);
  myBox.numParticles = myNumParticles;

  if(CkMyPe() == 0){
    CkPrintf("[STATS] node inter %lu part inter %lu open crit %lu dt %f\n", dtred->pnInteractions, dtred->ppInteractions, dtred->openCrit, globalParams.dtime);
  }

  iteration++;
  CkCallback cb;
  if(iteration == globalParams.iterations){
    cb = CkCallback(CkIndex_Main::niceExit(),mainProxy);
    contribute(0,0,CkReduction::sum_int,cb);
  }
  else if(iteration % globalParams.balancePeriod == 0){
    if(CkMyPe() == 0) CkPrintf("(%d) INITIATE LB\n", CkMyPe()); 
    for(int i = 0; i < submittedParticles.length()-1; i++){
      TreePiece *tp = submittedParticles[i].owner;
      tp->doAtSync();
    }
  }
  else{
    init();
    cb = CkCallback(CkIndex_DataManager::recvUnivBoundingBox(NULL),thisProxy);
    contribute(sizeof(BoundingBox),&myBox,BoundingBoxGrowReductionType,cb);
  }
  
  delete msg;
}

void DataManager::recvUnivBoundingBox(CkReductionMsg *msg){
  BoundingBox &univBB = *((BoundingBox *)msg->getData());
  decompose(univBB);
  delete msg;
}

void DataManager::freeCachedData(){
  map<Key,Request>::iterator it;

  for(it = particleRequestTable.begin(); it != particleRequestTable.end(); it++){
    Request &request = it->second;
    CkAssert(request.sent);
    CkAssert(request.data != NULL);
    CkAssert(request.requestors.length() == 0);
    CkAssert(request.msg != NULL);
    
    // no need to set the particles of a cached
    // bucket to NULL: we will delete the bucket
    // anyway
    if(!request.parentCached){
      request.parent->setParticles(NULL,0);
    }

    delete (ParticleReplyMsg *)(request.msg);
  }

  for(it = nodeRequestTable.begin(); it != nodeRequestTable.end(); it++){
    Request &request = it->second;
    CkAssert(request.sent);
    CkAssert(request.data != NULL);
    CkAssert(request.requestors.length() == 0);
    CkAssert(request.msg != NULL);

    // tell the root of the nodes in this 
    // fetched entry that its children don't
    // exist anymore
    // cached parents may be deleted before
    // their children, so we don't set their
    // children
    if(!request.parentCached){
      request.parent->setChildren(NULL,0);
    }

    delete (NodeReplyMsg *)(request.msg);
  }

  nodeRequestTable.clear();
  particleRequestTable.clear();
}

void DataManager::quiescence(){
  CkPrintf("QUIESCENCE dm %d pieces done %d (%d) nodereq %d partreq %d\n",
              CkMyPe(),
              numTreePiecesDoneTraversals,
              numLocalTreePieces,
              nodeReqs.test(),
              partReqs.test()
              );
  
  CkCallback cb(CkIndex_Main::quiescenceExit(),mainProxy);
  contribute(0,0,CkReduction::sum_int,cb);
}

void DataManager::freeTree(){
  if(root != NULL){
    FreeTreeWorker<ForceData> freeWorker;
    fillTrav.postorderTraversal(root,&freeWorker); 
    delete root;
    root = NULL;
  }
}

void DataManager::kickDriftKick(OrientedBox<Real> &box, Real &energy){
  Vector3D<Real> dv;

  Particle *pstart = myParticles.getVec();
  Particle *pend = myParticles.getVec()+myNumParticles;
  Real dt_k1, dt_k2;
  if(iteration == 0){
    dt_k1 = globalParams.dtime;
    // won't have KE stored by a previous 
    // iteration for this one
    for(Particle *p = pstart; p != pend; p++){
      energy += p->mass*(p->velocity.lengthSquared()); 
    }
    energy /= 2.0;
  }
  else{
    dt_k1 = globalParams.dthf;
    energy = savedEnergy;
    savedEnergy = 0.0;
  }

  dt_k2 = globalParams.dthf;

  for(Particle *p = pstart; p != pend; p++){
    energy += p->mass*p->potential;
    // kick
    p->velocity += dt_k1*p->acceleration;
    savedEnergy += p->mass*(p->velocity.lengthSquared()); 
    // drift
    p->position += globalParams.dtime*p->velocity;
    // kick
    p->velocity += dt_k2*p->acceleration;
    
    box.grow(p->position);

    p->acceleration = Vector3D<Real>(0.0);
    p->potential = 0.0;

  }
  savedEnergy /= 2.0;
}

void DataManager::findMinVByA(DtReductionStruct &dtred){
  if(myNumParticles == 0) {
    dtred.haveNaN = false;
    return;
  }
  
  dtred.haveNaN = false;

  for(int i = 0; i < myNumParticles; i++){
    Real v = myParticles[i].velocity.length();
    Real a = myParticles[i].acceleration.length();
    CkAssert(!isnan(v));
    if(isnan(a)) dtred.haveNaN = true;
  }

}

void DataManager::markNaNBuckets(){
#if 0
  for(int i = 0; i < myBuckets.length(); i++){
    Node<ForceData> *bucket = myBuckets[i];
    Particle *part = bucket->getParticles();
    int numParticles = bucket->getNumParticles();
    for(int j = 0; j < numParticles; j++){
      if(isnan(part[j].acceleration.length())){
        bucket->setType(Invalid);
        break;
      }
    }
  }
#endif
}

void DataManager::init(){
  CkAssert(pendingMoments.empty());
  // safe to reset here, since all tree pieces 
  // must have finished iteration
  freeCachedData();

  decompIterations = 0;

  haveRanges = false;
  numLocalTreePieces = -1;
  myBuckets.length() = 0;
  doneTreeBuild = false;
  treeMomentsReady = false;
  numTreePiecesDoneTraversals = 0;

  firstSplitterRound = true;
  freeTree();
  nodeTable.clear();
  localTPRoots.clear();

  CkAssert(activeBins.getNumCounts() == 0);

  if(CkMyPe() == 0 && keyRanges != NULL){
    delete[] keyRanges;
    keyRanges = NULL;
  }
  else if(CkMyPe() != 0){
    delete rangeMsg;
    rangeMsg = NULL;
  }


  numInteractions[0] = 0;
  numInteractions[1] = 0;
  numInteractions[2] = 0;
  submittedParticles.length() = 0;
}

void DataManager::resumeFromLB(){
  /* we delay the freeing of data structures to this point
   * because we want the tree pieces to be able to use the
   * tree, and in particular their roots for load balancing
   * */
  init();
  CkCallback cb(CkIndex_DataManager::recvUnivBoundingBox(NULL),thisProxy);
  contribute(sizeof(BoundingBox),&myBox,BoundingBoxGrowReductionType,cb);
}

#if 0
extern string NodeTypeColor[];
void DataManager::printTree(Node<ForceData> *nd, ostream &os){
  os << nd->getKey() 
     << "[label=\""<< nd->getKey() 
     << "," << nd->getNumParticles() 
     << "," << nd->getOwnerStart() << ":" << nd->getOwnerEnd() 
     << "\","
     << "style=\"filled\""
     << "color=\"" << NodeTypeColor[nd->getType()] << "\""
     << "]" << endl;
  if(nd->getOwnerStart() == nd->getOwnerEnd()) return;
  for(int i = 0; i < nd->getNumChildren(); i++){
    Node<ForceData> *child = nd->getChildren()+i;
    if(child != NULL){
      os << nd->getKey() << " -> " << child->getKey() << endl;
      printTree(child,os);
    }
  }
}

#endif

#include "Traversal_defs.h"

