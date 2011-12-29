#ifndef __WORKER_H__
#define __WORKER_H__

#include "Node.h"
#include "util.h"
#include "MultipoleMoments.h"
#include "Descriptor.h"

#include <map>
using namespace std;

// worker abstract class
template<typename T>
class CutoffWorker {
  public:
  virtual int work(Node<T> *) = 0; 
  virtual void work(ExternalParticle *) {}
  virtual void bucketDone(Key k) {}
  virtual void *getContext() {return NULL;}
  virtual void setContext(void *context) {}
  virtual void done() {}
  virtual void beforeParticleForces(Key k) {}
};

#if 0
class DataManager;
class ParticleFlushWorker : public CutoffWorker<NodeDescriptor> {
  int leafCnt;
  DataManager *dataManager;

  public:
  ParticleFlushWorker(DataManager *dm) : 
    dataManager(dm),
    leafCnt(0)
  {
  }

  int work(Node<NodeDescriptor> *node);
  int getNumLeaves(){
    return leafCnt;
  }
};

class MomentsWorker : public CutoffWorker<ForceData> {

  // tree pieces on this PE, in sorted order
  // this allows us to mark the nodes that are
  // internal to each PE. The last PE in this list
  // is a dummy value, and is equal to the number of
  // useful treepieces+1
  CkVec<TreePieceDescriptor> &peTreePieces;
  map<Key,Node<ForceData>*> &nodeTable;
  CkVec<Node<ForceData>*> &buckets;
  map<int,Node<ForceData>*> &tpRoots; 

  int curTP;

  public: 
  MomentsWorker(CkVec<TreePieceDescriptor> &petps,
                map<Key,Node<ForceData>*> &tab,
                map<int,Node<ForceData>*> &roots,
                CkVec<Node<ForceData>*> &bucks
                ) : 
    peTreePieces(petps),
    nodeTable(tab),
    tpRoots(roots),
    buckets(bucks),
    curTP(0)
  {
  }

  int work(Node<ForceData> *node);
  void setLeafType(Node<ForceData> *leaf);
  void setTypeFromChildren(Node<ForceData> *node);
};
#endif

class State;
class TraversalWorker : public CutoffWorker<ForceData> {
  protected:
  TreePiece *ownerTreePiece;
  State *state;
  
#ifdef SPLASH_COMPATIBLE
  Particle *currentParticle;
#else
  Node<ForceData> *currentBucket;
#endif

  TraversalWorker() : 
    ownerTreePiece(NULL),
#ifdef SPLASH_COMPATIBLE
    currentParticle(NULL)
#else
    currentBucket(NULL)
#endif
  {
  }

  public:

#ifdef SPLASH_COMPATIBLE
  void reset(TreePiece *owner, State *s, Particle *part)
#else
  void reset(TreePiece *owner, State *s, Node<ForceData> *bucket)
#endif
  {
    ownerTreePiece = owner;
    state = s;
#ifdef SPLASH_COMPATIBLE
    currentParticle = part;
#else
    currentBucket = bucket;
#endif
  }

  void *getContext(){
    void *ret;
#ifdef SPLASH_COMPATIBLE
    ret = (void *) currentParticle;
#else
    ret = (void *) currentBucket;
#endif
    return ret;
  }

  void setContext(void *context){
#ifdef SPLASH_COMPATIBLE
    currentParticle = (Particle *) context;
#else
    currentBucket = (Node<ForceData> *) context;
#endif
  }

  int work(Node<ForceData> *node);
  void work(ExternalParticle *particle);
  void bucketDone(Key k);
  
  virtual void done() {}
  virtual bool getKeep(NodeType type) = 0;
  virtual bool repeat(NodeType type) {return false;}

  void beforeParticleForces(Key k); 
};

class LocalTraversalWorker : public TraversalWorker {
  static const bool keep[];
  public:
  LocalTraversalWorker() : TraversalWorker() {}
  bool getKeep(NodeType type);
};

class RemoteTraversalWorker : public TraversalWorker {
  static const bool keep[];
  public:
  RemoteTraversalWorker() : TraversalWorker() {}
  void done();
  bool getKeep(NodeType type);
  bool repeat(NodeType type);
};

class TreeSizeWorker : public CutoffWorker<ForceData> {

  int numNodes;
  int cutoffDepth;

  public:
  TreeSizeWorker(int d) : 
    cutoffDepth(d),
    numNodes(0)
  {
  }

  int work(Node<ForceData> *node);
  int getNumNodes(){
    return numNodes;
  }
};

#endif
