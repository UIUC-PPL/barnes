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
  numTreePieces(1),
  firstSplitterRound(false),
  decompIterations(0),
  iteration(0),
  haveRanges(false),
  keyRanges(NULL),
  rangeMsg(NULL),
  numLocalTreePieces(-1),
  doneTreeBuild(false),
  treeMomentsReady(false),
  recvdAllUniqueNodes(false),
  numUniqueNodeMsgsRecvd(0),
  numTreePiecesDoneTraversals(0),
  prevIterationStart(0.0)
{
#ifdef STATISTICS
  numInteractions[0] = 0;
  numInteractions[1] = 0;
  numInteractions[2] = 0;
#endif
  savedEnergy = 0.0;
}

void DataManager::loadParticles(CkCallback &cb){
  numRankBits = LOG_BRANCH_FACTOR;

  const char *fname = globalParams.filename.c_str();
  int npart = globalParams.numParticles;


#ifdef MEMCHECK
  CkPrintf("[%d] before loadParticles\n", CkMyPe());
  CmiMemoryCheck();
#endif

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

  //CkPrintf("[%d] seeking to %d\n", myid, offset);
  //  CkAssert(partFile.is_open());
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

#ifdef DECOMP_VERY_VERBOSE
#endif

  partFile.close();

  contribute(sizeof(BoundingBox),&myBox,BoundingBoxGrowReductionType,cb);

#ifdef MEMCHECK
  CkPrintf("[%d] before loadParticles\n", CkMyPe());
  CmiMemoryCheck();
#endif
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

#ifdef CHECK_CORRUPTION
  CkPrintf("(%d) check before decompose\n", CkMyPe());
  CmiMemoryCheck();
#endif

  hashParticleCoordinates(universe.box);
  myParticles.quickSort();

  if(CkMyPe()==0){
    float memMB = (1.0*CmiMemoryUsage())/(1<<20);
    ostringstream oss; 
    /*
    CkPrintf("(%d) start iteration %d mem %f MB univ %f %f %f %f %f %f prev %g s\n", 
              CkMyPe(),iteration,memMB,
              universe.box.lesser_corner.x,
              universe.box.lesser_corner.y,
              universe.box.lesser_corner.z,
              universe.box.greater_corner.x,
              universe.box.greater_corner.y,
              universe.box.greater_corner.z,
              CmiWallTimer()-prevIterationStart
              );
    */
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
#ifdef MEMCHECK
  CkPrintf("[%d] before initHistogramParticles\n", CkMyPe());
  CmiMemoryCheck();
#endif

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

#ifdef MEMCHECK
  CkPrintf("[%d] after initHistogramParticles\n", CkMyPe());
  CmiMemoryCheck();
#endif
  //CkPrintf("[%d]Added %d particles\n", CkMyPe(), myNumParticles);
}

void DataManager::sendHistogram(){
  /*
  NodeDescriptor *descriptors = activeBins.getCounts();
  for(int i = 0; i < activeBins.getNumCounts(); i++){
    CkPrintf("(%d) sendHistogram nk %lu count %d sk %lx lk %lx\n",
              CkMyPe(), 
              descriptors->nodeKey,
              descriptors->numParticles,
              descriptors->smallestKey,
              descriptors->largestKey);
    CkAssert((descriptors->numParticles == 0) || (descriptors->smallestKey <= descriptors->largestKey));
    descriptors++;
  }
  */

#ifdef DECOMP_VERBOSE
#endif

  CkCallback cb(CkIndex_DataManager::receiveHistogram(NULL),0,this->thisgroup);
  contribute(sizeof(NodeDescriptor)*activeBins.getNumCounts(),activeBins.getCounts(),NodeDescriptorReductionType,cb);
  //contribute(sizeof(int)*activeBins.getNumCounts(),activeBins.getCounts(),CkReduction::sum_int,cb);
  activeBins.reset();

#ifdef MEMCHECK
  CkPrintf("[%d] after sendHistogram\n", CkMyPe());
  CmiMemoryCheck();
#endif
}

// executed on PE 0
void DataManager::receiveHistogram(CkReductionMsg *msg){

#ifdef MEMCHECK
  CkPrintf("[%d] before receiveHistogram\n", CkMyPe());
  CmiMemoryCheck();
#endif

  int numRecvdBins = msg->getSize()/sizeof(NodeDescriptor);
  NodeDescriptor *descriptors = (NodeDescriptor *)msg->getData();
#ifdef DECOMP_VERBOSE
  CkPrintf("PL0: in receiveHistogram %d bins\n", numRecvdBins);
#endif

  // XXX remove this and make a refine function for ActiveBinInfo
  CkVec<int> binsToRefine;

  binsToRefine.reserve(numRecvdBins);
  binsToRefine.length() = 0;

  int particlesHistogrammed = 0;

  CkVec<pair<Node<NodeDescriptor>*,bool> > *active = activeBins.getActive();
  CkAssert(numRecvdBins == active->length());

  for(int i = 0; i < numRecvdBins; i++){
#ifdef DECOMP_VERBOSE
    CkPrintf("PL0: bin %d nk %lu count %d sk %lx lk %lx\n", 
              i, descriptors[i].nodeKey, descriptors[i].numParticles, descriptors[i].smallestKey, descriptors[i].largestKey);
#endif
    if(descriptors[i].numParticles > (Real)(DECOMP_TOLERANCE*globalParams.ppc)){
      // need to refine this bin (partition)
      binsToRefine.push_back(i);
#ifdef DECOMP_VERBOSE
      CkPrintf("PL0: refine bin %d\n", i);
#endif
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

  int numBinsToRefine = binsToRefine.length();

  if(numBinsToRefine > 0){

#ifdef DECOMP_VERBOSE
    CkPrintf("DM0: Refine %d bins %d particles histogrammed\n", numBinsToRefine, particlesHistogrammed);
#endif
    SplitterMsg *m = new (numBinsToRefine,0) SplitterMsg;
    memcpy(m->splitBins,binsToRefine.getVec(),sizeof(int)*numBinsToRefine);
    m->nSplitBins = numBinsToRefine; 
    thisProxy.receiveSplitters(m);
    decompIterations++;
  }
  else{
    // create tree pieces and send proxy
#ifdef DECOMP_VERBOSE
#endif
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

#ifdef MEMCHECK
  CkPrintf("[%d] before loadParticles\n", CkMyPe());
  CmiMemoryCheck();
#endif

  delete msg;
}

void DataManager::flushParticles(){
  ParticleFlushWorker pfw(this);
  scaffoldTrav.preorderTraversal(sortingRoot,&pfw);

  int numUsefulTreePieces = pfw.getNumLeaves(); 
#ifdef DECOMP_VERBOSE
  if(CkMyPe() == 0) CkPrintf("(%d) num useful tree pieces %d\n", CkMyPe(), numUsefulTreePieces);
#endif

  for(int i = numUsefulTreePieces; i < globalParams.numTreePieces; i++){
    treePieceProxy[i].receiveParticles();
  }

  // done with sorting tree; delete
  FreeTreeWorker<NodeDescriptor> freeWorker;
  scaffoldTrav.postorderTraversal(sortingRoot,&freeWorker);

  delete sortingRoot;
  sortingRoot = NULL;

}

void DataManager::receiveSplitters(SplitterMsg *msg){

  int numRefineBins = msg->nSplitBins;
  /*
  if(firstSplitterRound){
    firstSplitterRound = false;
    freeTree();
    nodeTable.clear();
  }
  */

#ifdef DECOMP_VERBOSE
  CkPrintf("[%d] received %d splitters\n", CkMyPe(), numRefineBins);
#endif

#ifdef MEMCHECK
  CkPrintf("[%d] before loadParticles\n", CkMyPe());
  CmiMemoryCheck();
#endif

  // process bins to refine
  activeBins.processRefine(msg->splitBins,msg->nSplitBins);

  // We traverse the final tree to flush particles to 
  // appropriate tree pieces

  // here, we will know of bins that have not
  // been refined or deleted in the present 
  // iteration; we can send particles to these
  // CkVec<Node<NodeDescriptor>*> &unrefined = activeBins.getUnrefined();

#ifdef DECOMP_VERBOSE
  //CkPrintf("[%d] num UNREFINED %d\n", CkMyPe(), unrefined.length());
#endif

  //activeBins.processEmpty(msg->emptyBins,msg->nEmptyBins);

  sendHistogram();
  delete msg;

#ifdef MEMCHECK
  CkPrintf("[%d] before loadParticles\n", CkMyPe());
  CmiMemoryCheck();
#endif
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

#ifdef DECOMP_VERBOSE
  CkPrintf("[%d] send particles to %d tree pieces\n", CkMyPe(), numTreePieces);
#endif

  if(CkMyPe() != 0){
    numTreePieces = msg->numTreePieces;
    keyRanges = msg->keys;
    haveRanges = true;
    // delete this later
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


  // there are tree piece on this PE
  // these will eventually receive their respective
  // particles and submit them to the DM. Also, the
  // DM will receive the rangeKeys from the decomposition
  // leader (i.e. DM on PE 0) When all particles and 
  // rangeKeys have been received, the DM begins tree
  // construction

#ifdef PARTICLE_LOADER_VERBOSE
  CkPrintf("[%d] done with sendParticles\n", 
      CkMyPe());
#endif
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
  if(submittedParticles.length() == numLocalTreePieces &&
     haveRanges){
    processSubmittedParticles();
  }
}

void DataManager::processSubmittedParticles(){
  int offset = 0;

  submittedParticles.quickSort();
  
  myParticles.resize(myNumParticles);
  //CkPrintf("(%d) process %d submitted particles from %d local treepieces\n", CkMyPe(), myNumParticles, numLocalTreePieces);

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
    // discard node when:
    // 1. ownerEnd of node is < curTP
    // 2. if not 1, check for numparticles in node 
    // if ownerStart of node is > curTP, curTP = curTP->next
    // when a node its split, 
    // search_begin_idx = node->ownerStart*2;
    // search_end_idx = (node->ownerEnd*2+1)+1;
    // search_idx = binary_search_first_ge(Node::ParticleLevelKey(rightChildKey), tps, search_begin_idx, search_end_idx);
    // lchild->ownerStart = node->ownerStart;
    // lchild->ownerEnd = search_idx/2;
    // rchild->ownerStart = search_idx/2;
    // rchild->ownerEnd = node->ownerEnd;
    // if(EVEN(search_idx)) lchild->ownerEnd--;

    for(int i = 0; i < active->length(); i++){
      Node<ForceData> *node = (*active)[i].first;
      if((node->getOwnerEnd() > node->getOwnerStart()) || 
         (node->getNumParticles() > limit)){
        refines.push_back(i);
      }
    }

    abi.processRefine(refines.getVec(), refines.length());
    numFatNodes = abi.getNumCounts();
  }

}

void DataManager::makeMoments(){
  if(root == NULL) return;

  MomentsWorker mw(submittedParticles,
                   nodeTable,
                   uniqueNodes,
                   myBuckets
                   );
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
  for(int i = 0; i < replyTo.length(); i++){
    MomentsMsg *m = new (NUM_PRIORITY_BITS) MomentsMsg(node);
    *(int *)CkPriorityPtr(m) = RECV_MOMENTS_PRIORITY;
    CkSetQueueing(m,CK_QUEUEING_IFIFO);
    TB_DEBUG("(%d) responding to %d with node %lu\n", CkMyPe(), replyTo[i], node->getKey());
    thisProxy[replyTo[i]].receiveMoments(m);
  }
  replyTo.length() = 0;
}

void DataManager::receiveMoments(MomentsMsg *msg){
  Node<ForceData> *node = lookupNode(msg->data.key);
  CkAssert(node != NULL);

  // update moments of leaf and pass these on 
  // to parent recursively; if there are requests
  // for these nodes, respond to them
  updateLeafMoments(node,msg->data);
   
  delete msg;
}

void DataManager::exchangeUniqueNodes(){
  int num = uniqueNodes.length();
  MomentsExchangeMsg *msg = new (num,0) MomentsExchangeMsg;
  for(int i = 0; i < num; i++){
    Node<ForceData> *node = uniqueNodes[i];
    msg->data[i] = *node;
  }
  msg->numNodes = num;
  msg->fromPE = CkMyPe();
  // XXX make allgather/all-to-all bcast more efficient
  thisProxy.recvUniqueNodes(msg);
}

void DataManager::recvUniqueNodes(MomentsExchangeMsg *msg){
  if(!doneTreeBuild){
    uniqueNodeMsgs.push_back(msg); 
    return;
  }

  extendTree(msg);
  numUniqueNodeMsgsRecvd++;
  afterMomentsUpdated();
}

void DataManager::processMomentUpdates(){
  CkAssert(doneTreeBuild);
  for(int i = 0; i < uniqueNodeMsgs.length(); i++){
    MomentsExchangeMsg *msg = uniqueNodeMsgs[i];
    extendTree(msg);
  }
  numUniqueNodeMsgsRecvd += uniqueNodeMsgs.length();
  afterMomentsUpdated();
}

void DataManager::extendTree(MomentsExchangeMsg *msg){
  if(msg->fromPE != CkMyPe()){
    for(int i = 0; i < msg->numNodes; i++){
      MomentsExchangeStruct &data = msg->data[i];
      
      CkAssert(data.key > Key(0));

      Node<ForceData> *node = lookupNode(data.key);
      if(node == NULL){
        node = createNode(data.key);
        CkAssert(data.type != Boundary);
        if(data.type == Internal) node->setType(Remote);
        else if(data.type == Bucket) node->setType(RemoteBucket);
        else if(data.type == EmptyBucket) node->setType(RemoteEmptyBucket);
      }
      updateLeafMoments(node,data);
    }
  }
  delete msg;
}

Node<ForceData> *DataManager::createNode(Key k){
  Key curKey = Key(1);
  int curDepth = 0;

  if(root == NULL){
    root = new Node<ForceData>(curKey,curDepth,NULL,0);
    root->setType(Remote);
    nodeTable[curKey] = root;
  }

  Node<ForceData> *curNode = root;

  int kDepth = Node<ForceData>::getDepthFromKey(k);
  Key mask = Key((1<<LOG_BRANCH_FACTOR)-1);
  int nshift = ((kDepth-1)*LOG_BRANCH_FACTOR);
  mask <<= nshift;

  while(curKey != k){
    CkAssert(curDepth >= 0);
    CkAssert(curDepth < ((TREE_KEY_BITS-1)/LOG_BRANCH_FACTOR));

    if(curNode->getNumChildren() == 0){
      // need to create children of current node
      Node<ForceData> *child = new Node<ForceData>[BRANCH_FACTOR];
      curNode->setChildren(child,BRANCH_FACTOR);

      int childDepth = curDepth+1;
      Key childKey = (curKey << LOG_BRANCH_FACTOR);
      for(int i = 0; i < BRANCH_FACTOR; i++){
        child->setKey(childKey);
        child->setDepth(childDepth);
        child->setParent(curNode);
        child->setType(Remote);

        nodeTable[childKey] = child;
        childKey++;
        child++;
      }
    }

    // the children of this node (now) exist
    // they were created either during initial tree building
    // or by the code above

    int whichChild = ((k & mask)>>nshift);
    nshift -= LOG_BRANCH_FACTOR;
    mask >>= LOG_BRANCH_FACTOR;
    curNode = curNode->getChild(whichChild);
    curKey = ((curKey<<LOG_BRANCH_FACTOR)+whichChild);
    
    curDepth++;

  }

  return curNode;
}

void DataManager::updateLeafMoments(Node<ForceData> *node, MomentsExchangeStruct &data){
  //node->data.moments = moments;
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

void DataManager::afterMomentsUpdated(){
  if(numUniqueNodeMsgsRecvd == CkNumPes()){
    numUniqueNodeMsgsRecvd = 0;
    uniqueNodeMsgs.length() = 0;
    recvdAllUniqueNodes = true;
    if(treeMomentsReady){
      startTraversal();
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
    RequestMsg *msg = bufferedNodeRequests[i];
    requestNode(msg);
  }
  for(int i = 0; i < bufferedParticleRequests.length(); i++){
    RequestMsg *msg = bufferedParticleRequests[i];
    requestParticles(msg);
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
      descr.owner->prepare(root,myBuckets.getVec(),descr.bucketStartIdx,descr.bucketEndIdx);
      int tpIndex = descr.owner->getIndex();
      treePieceProxy[tpIndex].startTraversal();
      submittedParticles[i+1].bucketStartIdx = bucketIdx;
      start = bucketIdx;
    }
    TreePieceDescriptor &descr = submittedParticles[numLocalTreePieces-1];
    descr.bucketEndIdx = myBuckets.length();
    descr.owner->prepare(root,myBuckets.getVec(),descr.bucketStartIdx,descr.bucketEndIdx);
    int tpIndex = descr.owner->getIndex();
    treePieceProxy[tpIndex].startTraversal();
  }
  else if(numLocalTreePieces > 0){
    for(int i = 0; i < numLocalTreePieces; i++){
      TreePieceDescriptor &descr = submittedParticles[i];
      descr.owner->prepare(root,myBuckets.getVec(),0,0);
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
  if(!request.sent){
    partReqs.incrRequests();

    RequestMsg *reqMsg = new (NUM_PRIORITY_BITS) RequestMsg(key,CkMyPe());
    *(int *)CkPriorityPtr(reqMsg) = REQUEST_PARTICLES_PRIORITY;
    CkSetQueueing(reqMsg,CK_QUEUEING_IFIFO);

    CkAssert(leaf->getOwnerStart() == leaf->getOwnerEnd());
    int owner = leaf->getOwnerStart();
    treePieceProxy[owner].requestParticles(reqMsg);
    request.sent = true;
    
    request.parent = leaf;
    RRDEBUG("(%d) REQUEST particles %lu from tp %d\n", 
            CkMyPe(), key, owner);
  }
  request.requestors.push_back(Requestor(worker,state,traversal,worker->getContext()));
  partReqs.incrDeliveries();
}

void DataManager::requestParticles(RequestMsg *msg){
  if(!treeMomentsReady){
    bufferedParticleRequests.push_back(msg);
    return;
  }

  RRDEBUG("(%d) REPLY particles key %lu to %d\n", 
          CkMyPe(), msg->key, msg->replyTo);

  map<Key,Node<ForceData>*>::iterator it = nodeTable.find(msg->key);
  CkAssert(it != nodeTable.end());
  Node<ForceData> *bucket = it->second;
  CkAssert(bucket->getType() == Bucket);

  Particle *data = bucket->getParticles();
  int np = bucket->getNumParticles();

  ParticleReplyMsg *pmsg = new (np,NUM_PRIORITY_BITS) ParticleReplyMsg;
  *(int *)CkPriorityPtr(pmsg) = RECV_PARTICLES_PRIORITY;
  CkSetQueueing(pmsg,CK_QUEUEING_IFIFO);

  pmsg->key = msg->key;
  pmsg->np = np;
  for(int i = 0; i < np; i++){
    pmsg->data[i] = data[i];
  }

  thisProxy[msg->replyTo].recvParticles(pmsg);
  delete msg;
}

void DataManager::requestNode(Node<ForceData> *leaf, CutoffWorker<ForceData> *worker, State *state, Traversal<ForceData> *traversal){
  Key key = leaf->getKey();
  Request &request = nodeRequestTable[key];
  if(!request.sent){
    nodeReqs.incrRequests();

    RequestMsg *reqMsg = new (NUM_PRIORITY_BITS) RequestMsg(key,CkMyPe());
    *(int *)CkPriorityPtr(reqMsg) = REQUEST_NODE_PRIORITY;
    CkSetQueueing(reqMsg,CK_QUEUEING_IFIFO);

    int numOwners = leaf->getOwnerEnd()-leaf->getOwnerStart()+1;
    int requestOwner = leaf->getOwnerStart()+(rand()%numOwners);
    RRDEBUG("(%d) REQUEST node %lu from tp %d\n", 
            CkMyPe(), key, requestOwner);
    treePieceProxy[requestOwner].requestNode(reqMsg);
    request.sent = true;
    request.parent = leaf;
  }
  request.requestors.push_back(Requestor(worker,state,traversal,worker->getContext()));
  nodeReqs.incrDeliveries();
}

void DataManager::requestNode(RequestMsg *msg){
  if(!treeMomentsReady){
    bufferedNodeRequests.push_back(msg);
    return;
  }

  RRDEBUG("(%d) REPLY node %lu to %d\n", 
          CkMyPe(), msg->key, msg->replyTo);


  map<Key,Node<ForceData>*>::iterator it = nodeTable.find(msg->key);
  CkAssert(it != nodeTable.end());
  Node<ForceData> *node = it->second;

  if(node->getNumChildren() == 0){
    CkPrintf("[%d] children of leaf node %lu type %s requested!\n", CkMyPe(), node->getKey(), NodeTypeString[node->getType()].c_str());
    CkAbort("Leaf children request\n");
  }

  TreeSizeWorker tsz(node->getDepth()+globalParams.cacheLineSize);
  fillTrav.topDownTraversal_local(node,&tsz);

  int nn = tsz.getNumNodes();

  NodeReplyMsg *nmsg = new (nn,NUM_PRIORITY_BITS) NodeReplyMsg;
  *(int *)CkPriorityPtr(nmsg) = RECV_NODE_PRIORITY;
  CkSetQueueing(nmsg,CK_QUEUEING_IFIFO);

  nmsg->key = msg->key;
  nmsg->nn = nn;

  Node<ForceData> *emptyBuf = nmsg->data;
  node->serialize(NULL,emptyBuf,globalParams.cacheLineSize);
  CkAssert(emptyBuf == nmsg->data+nn);

  thisProxy[msg->replyTo].recvNode(nmsg);

  delete msg;
}

void DataManager::recvParticles(ParticleReplyMsg *msg){
  map<Key,Request>::iterator it = particleRequestTable.find(msg->key);
  CkAssert(it != particleRequestTable.end());

  RRDEBUG("(%d) RECVD particle REPLY for key %lu\n", 
          CkMyPe(), msg->key);

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

  RRDEBUG("(%d) RECVD node REPLY for key %lu\n", 
          CkMyPe(), msg->key);

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

  ostringstream oss;
  //ostream &oss = cerr;
  oss << "(" << CkMyPe() << ") key check: " << msg->key << endl;
  TreeChecker checker(oss);
  fillTrav.topDownTraversal_local(node,&checker);

  nodeReqs.decrRequests();
  nodeReqs.decrDeliveries(req.requestors.length());
  RRDEBUG("(%d) DELIVERING key %lu\n", 
          CkMyPe(), msg->key);

  req.deliverNode();
  RRDEBUG("(%d) DELIVERED key %lu\n", 
          CkMyPe(), msg->key);
}

#ifdef STATISTICS
void DataManager::traversalsDone(CmiUInt8 pnInter, CmiUInt8 ppInter, CmiUInt8 openCrit)
#else
void DataManager::traversalsDone()
#endif
{
  numTreePiecesDoneTraversals++;
#ifdef STATISTICS
  numInteractions[0] += pnInter;
  numInteractions[1] += ppInter;
  numInteractions[2] += openCrit;
#endif
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

  InteractionChecker ic;
  fillTrav.postorderTraversal(root,&ic);

  DtReductionStruct dtred;
  findMinVByA(dtred);

#ifdef STATISTICS
  dtred.pnInteractions = numInteractions[0];
  dtred.ppInteractions = numInteractions[1];
  dtred.openCrit = numInteractions[2];
  numInteractions[0] = 0;
  numInteractions[1] = 0;
  numInteractions[2] = 0;
#endif

  CkCallback cb(CkIndex_DataManager::advance(NULL),thisProxy);
  contribute(sizeof(DtReductionStruct),&dtred,DtReductionType,cb);

}

void DataManager::advance(CkReductionMsg *msg){

  DtReductionStruct *dtred = (DtReductionStruct *)(msg->getData());
  if(dtred->haveNaN){
    CkPrintf("(%d) iteration %d NaN accel detected! Exit...\n", CkMyPe(), iteration);
    markNaNBuckets();
    printTree();
    CkCallback exitCb(CkCallback::ckExit);
    contribute(0,0,CkReduction::sum_int,exitCb);
    return;
  }

  BoundingBox myBox;
  kickDriftKick(myBox.box,myBox.energy);

  Real pad = 0.001;
  myBox.expand(pad);
  myBox.numParticles = myNumParticles;

  if(CkMyPe() == 0){
#ifdef STATISTICS
    CkPrintf("[STATS] node inter %lu part inter %lu open crit %lu dt %f\n", dtred->pnInteractions, dtred->ppInteractions, dtred->openCrit, globalParams.dtime);
#endif
  }

  CkAssert(pendingMoments.empty());
  // safe to reset here, since all tree pieces 
  // must have finished iteration
  freeCachedData();

  submittedParticles.length() = 0;
  haveRanges = false;
  myBuckets.length() = 0;
  doneTreeBuild = false;
  numUniqueNodeMsgsRecvd = 0;
  uniqueNodes.length() = 0;
  recvdAllUniqueNodes = false;
  treeMomentsReady = false;
  numTreePiecesDoneTraversals = 0;

  firstSplitterRound = true;
  freeTree();
  nodeTable.clear();

  CkAssert(activeBins.getNumCounts() == 0);

  if(CkMyPe() == 0) delete[] keyRanges;
  else delete rangeMsg;

  iteration++;
  CkCallback cb;
  if(iteration == globalParams.iterations){
    cb = CkCallback(CkIndex_Main::niceExit(),mainProxy);
    contribute(0,0,CkReduction::sum_int,cb);
  }
  else{
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
  for(it = nodeRequestTable.begin(); it != nodeRequestTable.end(); it++){
    Request &request = it->second;
    CkAssert(request.sent);
    CkAssert(request.data != NULL);
    CkAssert(request.requestors.length() == 0);
    CkAssert(request.msg != NULL);

    // tell the root of the nodes in this 
    // fetched entry that its children don't
    // exist anymore
    request.parent->setChildren(NULL,0);

    delete (NodeReplyMsg *)(request.msg);
  }

  for(it = particleRequestTable.begin(); it != particleRequestTable.end(); it++){
    Request &request = it->second;
    CkAssert(request.sent);
    CkAssert(request.data != NULL);
    CkAssert(request.requestors.length() == 0);
    CkAssert(request.msg != NULL);

    request.parent->setParticles(NULL,0);

    delete (ParticleReplyMsg *)(request.msg);
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
  FreeTreeWorker<ForceData> freeWorker;
  fillTrav.postorderTraversal(root,&freeWorker); 
  if(root != NULL){
    delete root;
    root = NULL;
  }
}

void DataManager::printTree(){
  ostringstream oss;
  oss << "pe" << CkMyPe() << "." << iteration << ".dot";

  ofstream ofs(oss.str().c_str());

  ofs << "digraph PE" << CkMyPe() << "_" << iteration << " {" << endl;
  ofs << "node [style=\"filled\"]" << endl;
  if(root != NULL){
    Node<ForceData> &rootRef = *root;
    ofs << rootRef;
  }
  ofs << "}" << endl;
  ofs.close();
}

void DataManager::printBoundingBoxes(){
  ostringstream oss;
  oss << "bb.pe" << CkMyPe() << ".dat";
  ofstream ofs(oss.str().c_str());
  root->printBoundingBox(ofs);
  ofs.close();
}

void DataManager::addBucketNodeInteractions(Key k, CmiUInt8 pn){
#ifdef CHECK_NUM_INTERACTIONS
  Node<ForceData> *node = nodeTable[k];
  node->addNodeInteractions(pn);
#endif
}

void DataManager::addBucketPartInteractions(Key k, CmiUInt8 pp){
#ifdef CHECK_NUM_INTERACTIONS
  Node<ForceData> *node = nodeTable[k];
  node->addPartInteractions(pp);
#endif
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

void DataManager::resetParticleAccelerations(){
  for(int i = 0; i < myNumParticles; i++){
    myParticles[i].acceleration.x = 0.0;
    myParticles[i].acceleration.y = 0.0;
    myParticles[i].acceleration.z = 0.0;
  }
}

void DataManager::findMinVByA(DtReductionStruct &dtred){
  if(myNumParticles == 0) {
    dtred.haveNaN = false;
    dtred.vbya = -1.0; 
    return;
  }
  
  /*
  Real min_v = myParticles[0].velocity.length();
  Real min_a = myParticles[0].acceleration.length();
  */
  dtred.haveNaN = false;

  //CkAssert(!isnan(min_v));
  //if(isnan(min_a)) dtred.haveNaN = true;

  //Real min_vbya = min_v/min_a;

  for(int i = 0; i < myNumParticles; i++){
    Real v = myParticles[i].velocity.length();
    Real a = myParticles[i].acceleration.length();
    CkAssert(!isnan(v));
    /*
    if(isnan(a)){
      CkPrintf("(%d) particle %d key %lu accel is nan!\n", CkMyPe(), i, myParticles[i].key);
      CkAbort("bad acceleration\n");
    }
    */
    if(isnan(a)) dtred.haveNaN = true;
    /*
    Real vbya = v/a;
    if(vbya < min_vbya){
      min_vbya = vbya;
      min_v = v; 
      min_a = a; 
    }
    */
  }

  //CkPrintf("(%d) vbya %f v %f a %f\n", CkMyPe(), min_vbya, min_v, min_a);
  //dtred.vbya = min_vbya;
  dtred.vbya = -1.0;
}

void DataManager::markNaNBuckets(){
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
}



#include "Traversal_defs.h"

