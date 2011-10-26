#ifndef __DATA_MANAGER_H__
#define __DATA_MANAGER_H__

#include "Particle.h"

#include "OrientedBox.h"
#include "barnes.decl.h"
#include "Node.h"
#include "Descriptor.h"
#include "ActiveBinInfo.h"

#include "Traversal_decls.h"
#include "Request.h"

class TreePiece;

#include <map>
using namespace std;

class TreePieceCounter : public CkLocIterator {            
  public:
  int count;
  CkHashtableT<CkArrayIndex, int> registered;               
  TreePieceCounter() : count(0) { }                      
  void addLocation(CkLocation &loc) {
    registered.put(loc.getIndex()) = ++count;               
  }
  void reset() {                                            
    count = 0;
    registered.empty();                                     
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

class DataManager : public CBase_DataManager {
  int numRankBits;
  double prevIterationStart;

  CkVec<Particle> myParticles;
  int myNumParticles;
  BoundingBox myBox;

  bool firstSplitterRound;

  int numTreePieces;
  int numLocalUsefulTreePieces;
  int numLocalTreePieces;
  TreePieceCounter localTreePieces;
  CkVec<TreePieceDescriptor> submittedParticles;

  int iteration;
  int decompIterations;
  ActiveBinInfo<ForceData> activeBins;
  // Root of particles on this PE
  // Used in both decomposition and traversal
  Node<ForceData> *root;


  CkVec<Node<ForceData>*> myBuckets;

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
  CkVec< std::pair<Key, int> > bufferedNodeRequests;
  CkVec< std::pair<Key, int> > bufferedParticleRequests;

  map<Key,Request> nodeRequestTable;
  map<Key,Request> particleRequestTable;

  Traversal<ForceData> fillTrav;
  int numTreePiecesDoneTraversals;
  CacheStats nodeReqs;
  CacheStats partReqs;

  /* used only on PE 0 to track drift in energy */ 
  Real compareEnergy;

  CmiUInt8 numInteractions[3];

  void kickDriftKick(OrientedBox<Real> &box, Real &energy);

  void hashParticleCoordinates(OrientedBox<Real> &universe);
  void initHistogramParticles();
  void sendHistogram();
  
  void senseTreePieces();
  void processSubmittedParticles();
  int buildTree(Node<ForceData> *node, int pstart, int pend, int tpstart, int tpend);
  void singleBuildTree(Node<ForceData> *node, TreePieceDescriptor &ownerTreePiece);

  void flushParticles();
  int flushAndMark(Node<ForceData> *node, int mark); 

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

  void startTraversal();
  void flushBufferedRemoteDataRequests();

  void freeCachedData();
  void freeTree();
  void finishIteration();

#if 0
  void findMinVByA(DtReductionStruct &);
  void markNaNBuckets();
#endif
  
  void printTree(Node<ForceData>*, ostream &);
  void doPrintTree();

  void init();

  public:
  DataManager();

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

  // called by tree piece that is making a request
  void requestNode(Node<ForceData> *leaf, CutoffWorker<ForceData> *worker, State *state, Traversal<ForceData> *callbackTraversal);
  void requestParticles(Node<ForceData> *leaf, CutoffWorker<ForceData> *worker, State *state, Traversal<ForceData> *callbackTraversal);

  // called by tree piece that is forwarding a remote request
  void requestNode(std::pair<Key, int> request);
  void requestParticles(std::pair<Key, int> request);
  
  void recvParticles(ParticleReplyMsg *msg);
  void recvNode(NodeReplyMsg *msg);

  void recvUnivBoundingBox(CkReductionMsg *msg);

  void quiescence();

  void addBucketNodeInteractions(Key k, CmiUInt8 pn);
  void addBucketPartInteractions(Key k, CmiUInt8 pp);

  void resumeFromLB();
};

#endif
