#ifndef __STATE_H__
#define __STATE_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>

using namespace std;

#if 0
extern string NodeTypeString[];
#endif

class TreePiece;
struct State {
  int pending;
  int current;
  Node<ForceData> **currentBucketPtr;
  TreePiece *ownerTreePiece;

#ifdef DEBUG_TRAVERSALS
  map<Key,set<Key> > bucketKeys;
#endif

  CmiUInt8 numInteractions[3];

  void incrPending() { pending++; }

  void decrPending(int n=1){
    pending -= n;
  }

  bool complete() { return (pending == 0); }

  State() : 
    pending(-1), 
    current(-1), 
    currentBucketPtr(NULL),
    ownerTreePiece(NULL)
  {

    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;
  }
  
  void reset(TreePiece *owner, int p, Node<ForceData> **bucketPtr){
  //void reset(TreePiece *owner, string &id, int p, Node<ForceData> **bucketPtr){
    ownerTreePiece = owner;
    //description = id;
    pending = p;
    current = 0;
    currentBucketPtr = bucketPtr;
  }

  void nodeEncountered(Key bucketKey, Node<ForceData> *node);
  void nodeOpened(Key bucketKey, Node<ForceData> *node);
  void nodeDiscarded(Key bucketKey, Node<ForceData> *node);
  void nodeComputed(Node<ForceData> *bucket, Key nodeKey);
  void bucketComputed(Node<ForceData> *bucket, Key k);

  void finishedIteration(){
    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;
  }

  void insert(Key bucketKey, Key k){
#ifdef DEBUG_TRAVERSALS
    Key parentKey = (k>>1);
    set<Key>::iterator it = bucketKeys[bucketKey].find(parentKey);
    while(it != bucketKeys[bucketKey].end()){
      bucketKeys[bucketKey].erase(it);
      parentKey >>= 1;
      it = bucketKeys[bucketKey].find(parentKey);
    }
    bucketKeys[bucketKey].insert(parentKey);
#endif
  }

  void incrPartNodeInteractions(Key bucketKey, CmiUInt8 n){
    numInteractions[0] += n;
  }

  void incrPartPartInteractions(Key bucketKey, CmiUInt8 n){
    numInteractions[1] += n;
  }

  void incrOpenCriterion(){
    numInteractions[2]++;
  }

  virtual string getDescription() = 0;

#if 0
  CkVec<Vector3D<Real> > savedAccelerations;
#endif
  void beforeForces(Node<ForceData> *bucket, Key k);

};

struct LocalState : public State {
  string getDescription(){ return "LOCAL"; }
};

struct RemoteState : public State {
  string getDescription(){ return "REMOTE"; }
};

#endif
