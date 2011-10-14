#ifndef __STATE_H__
#define __STATE_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <map>

using namespace std;

extern string NodeTypeString[];

//typedef void (*StateCompletionFn)(void *context);

class TreePiece;
struct State {
  int pending;
  int current;
  Node<ForceData> **currentBucketPtr;
  TreePiece *ownerTreePiece;

#ifdef DEBUG_TRAVERSALS
  map<Key,set<Key> > bucketKeys;
#endif

#ifdef STATISTICS
  CmiUInt8 numInteractions[3];
#endif

#ifdef TREE_PIECE_LOG
  ofstream logFile;
#endif

#ifdef CHECK_NUM_INTERACTIONS
  map<Key,CmiUInt8> bucketNodeInteractions;
  map<Key,CmiUInt8> bucketPartInteractions;
#endif

  string description;

  void incrPending() { pending++; }
  bool decrPending() {
    pending--;
    return complete();
  }

  bool decrPending(int n){
    pending -= n;
    return complete();
  }

  bool complete() { return (pending == 0); }

  State() : 
    pending(-1), 
    current(-1), 
    currentBucketPtr(NULL),
    ownerTreePiece(NULL)
  {

#ifdef STATISTICS
    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;
#endif
  }
  
  State(int i, string id, TreePiece *tp) : 
    pending(-1), 
    current(-1), 
    ownerTreePiece(tp),
    currentBucketPtr(NULL) 
  {
#ifdef TREE_PIECE_LOG
    ostringstream oss;
    oss << "tp." << i << "." << id << ".log";
    logFile.open(oss.str().c_str());
#endif
    description = id;

#ifdef STATISTICS
    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;
#endif
  }

  void reset(int p, Node<ForceData> **bucketPtr){
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
#ifdef TREE_PIECE_LOG
    logFile.close();
#endif
#ifdef STATISTICS
    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;
#endif
#ifdef CHECK_NUM_INTERACTIONS
    bucketNodeInteractions.clear();
    bucketPartInteractions.clear();
#endif
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
#ifdef STATISTICS
    numInteractions[0] += n;
#endif
#ifdef CHECK_NUM_INTERACTIONS
    bucketNodeInteractions[bucketKey] += n;
#endif
  }

  void incrPartPartInteractions(Key bucketKey, CmiUInt8 n){
#ifdef STATISTICS
    numInteractions[1] += n;
#endif
#ifdef CHECK_NUM_INTERACTIONS
    bucketPartInteractions[bucketKey] += n;
#endif
  }

  void incrOpenCriterion(){
#ifdef STATISTICS
    numInteractions[2]++;
#endif
  }
};

/*
struct TraversalState {
  int curBucket;

  TraversalState() : curBucket(0), State(cb,context,pending) {}
};
*/


#endif
