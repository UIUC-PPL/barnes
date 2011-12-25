#ifndef __DATA_MANAGER_H__
#define __DATA_MANAGER_H__

#include "Particle.h"

#include "OrientedBox.h"
#include "MeshStreamer.h"
#include "barnes.decl.h"
#include "Node.h"
#include "Descriptor.h"
#include "ActiveBinInfo.h"

#include "Traversal_decls.h"
#include "Request.h"


class TreePiece;
extern CProxy_TreePiece treePieceProxy;

#include <map>
using namespace std;

class DataManager;
class TreePieceCounter : public CkLocIterator {            
  public:
  int count;
  int numParticles;
  CkVec<TreePieceDescriptor> submittedParticles;

  TreePieceCounter() { }                      

  void addLocation(CkLocation &loc);
  void reset() {                                            
    numParticles = 0;
    count = 0;
    submittedParticles.length() = 0;
  }                                                         
};

struct CacheStats {
  int outstandingRequests;
  int outstandingDeliveries;

  CacheStats() : 
    outstandingRequests(0),
    outstandingDeliveries(0)
  {
  }

  void incrRequests(){ outstandingRequests++; }
  void decrRequests(int n=1){ outstandingRequests -= n; }
  void incrDeliveries(){ outstandingDeliveries++; }
  void decrDeliveries(int n=1){ outstandingDeliveries -= n; }
  bool test(){ return (outstandingRequests==0) && (outstandingDeliveries==0); }
};

class DataManager : public MeshStreamerClient<NodeRequest> {
  int numRankBits;
  int numPesPerNode;
  double prevIterationStart;

  double phaseTime;

  CkVec<Particle> myParticles;
  int myNumParticles;
  BoundingBox myBox;

  int numTreePieces;
  int numLocalUsefulTreePieces;
#if 0
  int numLocalTreePieces;
#endif
  TreePieceCounter localTreePieces;
#if 0
  //CkVec<TreePieceDescriptor> submittedParticles;
#endif

  int numMomentsRequested;
  int numMomentsReceived;

  int iteration;
  int decompIterations;
  ActiveBinInfo<ForceData> activeBins;
  // Root of particles on this PE
  // Used in both decomposition and traversal
  Node<ForceData> *root;


  CkVec<Node<ForceData>*> myBuckets;

#if 0
  bool doneFlushParticles; 
#endif
  // I am done constructing the tree 
  // from particles present on this PE
  bool doneTreeBuild;
  map<Key,Node<ForceData>*> nodeTable;

  map<Key,CkVec<int> > pendingMoments;

  // I have processed the moment 
  // contributions from all other PEs, so that
  // the tree on this PE is now ready for 
  // traversal
  bool treeMomentsReady;
  //CkVec<RequestMsg*> bufferedNodeRequests;
  CkVec< std::pair<Key, int> > bufferedNodeRequests;
  CkVec< std::pair<Key, int> > bufferedParticleRequests;

  map<Key,Request> nodeRequestTable;
  map<Key,Request> particleRequestTable;

  // needed for skipping decomposition
  CkVec<Key> treePieceRoots;
  bool doSkipDecomposition;

  Traversal<ForceData> fillTrav;
  int numTreePiecesDoneTraversals;
  int numTreePiecesDoneRemoteRequests;
  CacheStats nodeReqs;
  CacheStats partReqs;

  /* used only on PE 0 to track drift in energy */ 
  Real compareEnergy;

  CmiUInt8 numInteractions[3];
  CProxy_DataManager myProxy;
  
  MeshStreamer<NodeRequest> *combiner;
  CkArray *tpArray;

  // Called by requestNode() and process() (when streaming-combining)
  void processNodeRequest(Key key, int replyTo);

  void kickDriftKick(OrientedBox<double> &box, Real &energy);

  void hashParticleCoordinates(OrientedBox<double> &universe);
  void initHistogramParticles();
  void sendHistogram();
  
  void senseTreePieces();
  void buildTree(Node<ForceData> *node, int pstart, int pend, int tpstart, int tpend);
  void singleBuildTree(Node<ForceData> *node, TreePieceDescriptor &ownerTreePiece);

  void flushParticles();
  int flushAndMark(Node<ForceData> *node, int mark); 
  // for skipping decomposition
  void skipFlushParticles();

  // helper function
  void sendParticleMsg(int tp, Particle *p, int numParticles);

#if 0
  void makeMoments();
  void flushMomentRequests();
#endif

  void notifyParentMomentsDone(Node<ForceData> *node);
  void momentsResponseHelper(Node<ForceData> *);
  void respondToMomentsRequest(Node<ForceData> *,CkVec<int>&);
  Node<ForceData> *lookupNode(Key k);

  void updateLeafMoments(Node<ForceData> *node, MomentsExchangeStruct &data);
  void passMomentsUpward(Node<ForceData> *node);
  void treeReady();

  void flushBufferedRemoteDataRequests();

  void freeCachedData();
  // to free the local tree at the end of the iteration, when not skipping decomposition
  void freeTree();
  // if skipping decomposition, we can reuse tree that already exists: reset its data 
  // (particles,type,moments,numChildrenMomentsReady
  void reuseTree();
  void finishIteration();

#if 0
  void findMinVByA(DtReductionStruct &);
  void markNaNBuckets();
#endif
  
  void printTree(Node<ForceData>*, ostream &);
  void doPrintTree(string name);

  void init();

  void nodeLevelMerge();
  void registerTopLevelNodes(Node<ForceData> *node, int tpstart, int tpend);

  public:
  DataManager();
  DataManager(CkMigrateMessage *m) {}

  void loadParticles(CkCallback &cb);

  void decompose(BoundingBox &universe);
  void receiveHistogram(CkReductionMsg *msg);
  void receiveSplitters(CkVec<int> splitBins);
  void sendParticles(int);
  void sendParticlesToTreePiece(Node<ForceData> *nd, int tp);

  void receiveMoments(MomentsExchangeStruct moments);
  
  // called by tree pieces
  void submitParticles(CkVec<ParticleMsg *> *vec, int numParticles, TreePiece *tp);
  void requestMoments(Key k, int replyTo);
  void advance(CkReductionMsg *);
  void traversalsDone(CmiUInt8 pnInter, CmiUInt8 ppInter, CmiUInt8 openCrit);
  void doneRemoteRequests();

  // called by tree piece that is making a request
  Node<ForceData> *requestNode(Node<ForceData> *leaf, CutoffWorker<ForceData> *worker, State *state, Traversal<ForceData> *callbackTraversal);
  ExternalParticle * requestParticles(Node<ForceData> *leaf, CutoffWorker<ForceData> *worker, State *state, Traversal<ForceData> *callbackTraversal);

  void combineNodeRequest(int tpindex, Key k);
  void combineParticleRequest(Key k, int tpindex);

  // after detecting quiescence, identify tree pieces on PE
  // and get their particles
  void processSubmittedParticles();

  void startTraversal();
  void doneNodeLevelMerge(PointerContainer);

  // called by tree piece that is forwarding a remote request
  //void requestNode(RequestMsg *);
  void requestNode(std::pair<Key, int> &request);
  void requestParticles(std::pair<Key, int> &request);
  
  void recvParticles(ParticleReplyMsg *msg);
  void recvNode(NodeReplyMsg *msg);

  void recvUnivBoundingBox(CkReductionMsg *msg);

  void quiescence();

  void addBucketNodeInteractions(Key k, CmiUInt8 pn);
  void addBucketPartInteractions(Key k, CmiUInt8 pp);

  void resumeFromLB();

  void process(NodeRequest &);
  void pup(PUP::er &p);

  // phase timing measurements
  void doneParticleFlush();
  void resumeProcessSubmittedParticles();
  void doneTreeBuilds();
  void resumeDoneNodeLevelMerge();
  void doneForces();
  void resumeFinishIteration();
};

#endif
