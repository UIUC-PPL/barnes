#ifndef __TREE_PIECE_H__
#define __TREE_PIECE_H__

#include "defines.h"
#include "OrientedBox.h"
#include "NDMeshStreamer.h"
#include "barnes.decl.h"
#include "Particle.h"
#include "Messages.h"
#include "Node.h"
#include "Vector3D.h"
#include "State.h"

#include "Worker.h"

#include "MultipoleMoments.h"
#include "Traversal_decls.h"

extern CProxy_DataManager dataManagerProxy;

class TreePiece : public MeshStreamerArrayClient<NodeRequest> {

  CProxy_TreePiece myProxy;
  CkVec<ParticleMsg *> decompMsgsRecvd;
  int myNumParticles;

  int iteration;

#ifndef SPLASH_COMPATIBLE
  int myNumBuckets;
  Node<ForceData> **myBuckets;
#endif
  Node<ForceData> *root;
  Node<ForceData> *myRoot;

#ifdef DEBUG_TRAVERSALS
  map<Key,set<Key> > bucketKeys;
#endif

  LocalState localTraversalState;
  RemoteState remoteTraversalState;

  LocalTraversalWorker localTraversalWorker;
  RemoteTraversalWorker remoteTraversalWorker;

  Traversal<ForceData> trav;

  DataManager *myDM;

  CkGroupID orbLBProxy;
  int numLB;
  bool haveOrbLB;

  int totalNumTraversals;

  void finishIteration();
  void init();
  void findOrbLB();

  void checkTraversals();
  void clearBucketsDebug();

  public:
  TreePiece();
  TreePiece(CkMigrateMessage *) {}

  inline int getIndex() {return thisIndex.data[0];}

  void receiveParticles(ParticleMsg *msg);

  // invoked by DM to get particles
  int getNumParticles() {return myNumParticles;}
  CkVec<ParticleMsg*> *getBufferedParticleMsgs() {return &decompMsgsRecvd;}

#ifdef SPLASH_COMPATIBLE
  void prepare(Node<ForceData> *_root, Node<ForceData> *_myRoot);
#else
  void prepare(Node<ForceData> *_root, Node<ForceData> *_myRoot, Node<ForceData> **buckets, int numBuckets);
#endif

  void startTraversal(int dummy);

  void doLocalGravity(RescheduleMsg *);
  void doRemoteGravity(RescheduleMsg *);

  void requestParticles(Key k, int replyTo);
  //void requestNode(RequestMsg *);
  void requestNode(Key k, int replyTo);

  void requestMoments(Key k, int replyTo);
  void traversalDone();
  // Promise the streaming library that we won't use it anymore
  void doneRemoteRequests();

  void quiescence();
  int getIteration();

  void pup(PUP::er &p);
  void startlb();
  void cleanup();
  void doAtSync(); 
  void ResumeFromSync();

  void process(const NodeRequest &req);
};

#endif
