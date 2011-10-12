#ifndef __DESCRIPTOR_H__
#define __DESCRIPTOR_H__

#include "defines.h"
#include "OrientedBox.h"
#include "MultipoleMoments.h"
#include "charm++.h"

struct BoundingBox {
  OrientedBox<Real> box;
  int numParticles;

  BoundingBox(){
    reset();
  }

  void reset(){
    numParticles = 0;
    box.reset();
  }

  void grow(const Vector3D<Real> &v){
    box.grow(v);
  }

  void grow(const BoundingBox &other){
    if(other.numParticles == 0) return;
    if(numParticles == 0){
      *this = other;
    }
    else{
      box.grow(other.box);
      numParticles += other.numParticles;
    }
  }

  void expand(Real pad){
    box.greater_corner = box.greater_corner*pad+box.greater_corner;
    box.lesser_corner = box.lesser_corner-pad*box.lesser_corner;
  }

  void pup(PUP::er &p){
    p | box;
    p | numParticles;
  }

};

// how many particles do I have under this node?
// what are the largest and smallest keys among them?
struct NodeDescriptor {
  int numParticles;
  Key nodeKey;

  Key smallestKey;
  Key largestKey;

  NodeDescriptor() {
    numParticles = 0;
    largestKey = Key(0);
    smallestKey = ~largestKey;
    nodeKey = smallestKey;
  }

  NodeDescriptor(int np, Key nk, Key sk, Key lk) : 
    numParticles(np), 
    nodeKey(nk), smallestKey(sk), largestKey(lk) 
  {
  }

  /*
  NodeDescriptor(Key nk) :
    numParticles(0), 
    nodeKey(nk) 
  {
    largestKey = Key(0);
    smallestKey = ~largestKey;
  }
  */

  void grow(const NodeDescriptor &other){
    /*
    CkPrintf("(%d) grow1 %d sk %lx lk %lx grow2 %d sk %lx lk %lx\n",
              CkMyPe(),
              numParticles, 
              smallestKey, 
              largestKey,
              other.numParticles, 
              other.smallestKey, 
              other.largestKey
              );
    */
    if(other.numParticles == 0){
      //CkPrintf("(%d) growresult 1\n", CkMyPe());
      return;
    }

    if(numParticles == 0){
      smallestKey = other.smallestKey;
      largestKey = other.largestKey;
      numParticles = other.numParticles;
      //CkPrintf("(%d) growresult 2\n", CkMyPe());
    } else {
      if(smallestKey > other.smallestKey) smallestKey = other.smallestKey;
      if(largestKey < other.largestKey) largestKey = other.largestKey;
      numParticles += other.numParticles;
      //CkPrintf("(%d) growresult 3\n", CkMyPe());
    }

  }
};

class TreePiece;
struct ParticleMsg;
struct TreePieceDescriptor {
  CkVec<ParticleMsg*> *vec;
  TreePiece *owner;
  int index;
  int numParticles;
  Key smallestKey;
  Key largestKey;

  int bucketStartIdx;
  int bucketEndIdx;

  TreePieceDescriptor() : 
    vec(NULL), owner(NULL), numParticles(0), index(-1), smallestKey(~Key(0)), largestKey(Key(0))
  {
  }

  TreePieceDescriptor(CkVec<ParticleMsg*> *v, int np, TreePiece *o, int i, Key sk, Key lk) : 
    vec(v), owner(o), index(i), smallestKey(sk), largestKey(lk), numParticles(np)
  {
  }

  TreePieceDescriptor(int i) : 
    vec(NULL), owner(NULL), numParticles(0), index(i), smallestKey(~Key(0)), largestKey(Key(0))
  {
  }

  bool operator<=(const TreePieceDescriptor &t){
    return smallestKey <= t.smallestKey;
  }
  bool operator>=(const TreePieceDescriptor &t){
    return smallestKey >= t.smallestKey;
  }
};

struct ForceData {
  OrientedBox<Real> box;
  MultipoleMoments moments;

  ForceData() 
  {
  }
};

struct DtReductionStruct {
#ifdef STATISTICS 
  CmiUInt8 pnInteractions;
  CmiUInt8 ppInteractions;
#endif
  Real vbya;
  bool haveNaN;

  DtReductionStruct &operator+=(const DtReductionStruct &other){
#ifdef STATISTICS
    pnInteractions += other.pnInteractions;
    ppInteractions += other.ppInteractions;
#endif
    if(other.haveNaN) haveNaN = true;

    if(other.vbya < 0.0) {} 
    else if(vbya < 0.0) vbya = other.vbya;
    else if(other.vbya < vbya) vbya = other.vbya;

    return *this;
  }
};


#endif
