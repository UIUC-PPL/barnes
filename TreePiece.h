#ifndef __TREE_PIECE_H__
#define __TREE_PIECE_H__

#include "defines.h"
#include "OrientedBox.h"
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

class TreePiece : public CBase_TreePiece {
  int numDecompMsgsRecvd;
  CkVec<ParticleMsg *> decompMsgsRecvd;
  int myNumParticles;

  int iteration;

  int myNumBuckets;
  Node<ForceData> **myBuckets;
  Node<ForceData> *root;
  Node<ForceData> *myRoot;

  LocalState localTraversalState;
  RemoteState remoteTraversalState;

  LocalTraversalWorker localTraversalWorker;
  RemoteTraversalWorker remoteTraversalWorker;

  Traversal<ForceData> trav;

  Key smallestKey;
  Key largestKey;
  DataManager *myDM;

  CkGroupID orbLBProxy;
  int numLB;
  bool haveOrbLB;

  int totalNumTraversals;

  void finishIteration();
  void init();
  void findOrbLB();

  public:
  TreePiece();
  TreePiece(CkMigrateMessage *) {}

  int getIndex() {return thisIndex;}

  void receiveParticles(ParticleMsg *msg);
  //void receiveParticles();
  void submitParticles();

  void prepare(Node<ForceData> *_root, Node<ForceData> *_myRoot, Node<ForceData> **buckets, int bucketStart, int bucketEnd);
  void startTraversal();

  void doLocalGravity(RescheduleMsg *);
  void doRemoteGravity(RescheduleMsg *);

  void requestParticles(std::pair<Key, int> request);
  void requestNode(std::pair<Key, int> request);

  void requestMoments(Key k, int replyTo);
  void traversalDone();

  void quiescence();
  int getIteration();

  void pup(PUP::er &p);
  void startlb();
  void doAtSync(); 
  void ResumeFromSync();
};

#endif
