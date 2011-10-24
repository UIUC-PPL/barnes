#ifndef __STATE_H__
#define __STATE_H__

/*
 * CharmBH: State.h
 * Track the state of a traversal. This structure records the
 * number of outstanding requests, and whether all buckets 
 * have been processed for a particular type of traversal.
 * Each traversal has a separate State.
 */

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
  /*
   * Num. pending buckets; initially set to number of 
   * buckets owned by TreePiece, to indicate that traversals
   * need to be initiated for those buckets. Incremented by Traversal
   * whenever it encounters a remote node whose children
   * must be fetched, or a remote bucket whose particles must
   * be fetched. When nodes/particles are received, this
   * count is decreased. Also decremented when a TreePiece 
   * finishes a traversal with a bucket.
   */
  int pending;
  // Index of current bucket with which tree is being traversed
  // This field is used in TreePiece::doLocalGravity and 
  // TreePiece::doRemoteGravity 
  int current;
  // Pointer to current bucket 
  Node<ForceData> **currentBucketPtr;
  // The tree piece that initiated the traversal being tracked
  TreePiece *ownerTreePiece;

#ifdef DEBUG_TRAVERSALS
  map<Key,set<Key> > bucketKeys;
#endif

  // Number of interactions and opening tests performed during the traversal
  CmiUInt8 numInteractions[3];

  // Called when the traversal misses on remote data
  void incrPending() { pending++; }

  // Either a remote request has been fulfilled or a bucket has finished traversal
  void decrPending(int n=1){
    pending -= n;
  }

  // A traversal is complete when it has no pending remote data requests and
  // all buckets have finished traversal
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
  
  // Reinitialize state
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

  // Clear counts of interactions/opening criterion applications
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

  // Methods to increment different statistics
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

#ifdef VERBOSE_TRAVERSAL_INTERACTION
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
