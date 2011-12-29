#ifndef __DESCRIPTOR_H__
#define __DESCRIPTOR_H__

/*
 * CharmBH - Descriptor.h
 * Descriptor structures that are used in various parts of the code.
 */


#include "defines.h"
#include "OrientedBox.h"
#include "MultipoleMoments.h"
#include "charm++.h"

/*
 * BoundingBox:
 * Used to calculate the bounding box of all particles in the
 * simulation universe. It also keeps track of particle energy,
 * to see whether there is a drift in the total system energy 
 * through the simulation.
 */
struct BoundingBox {
  OrientedBox<double> box;
  int numParticles;
  Real energy;

  BoundingBox(){
    reset();
  }

  void reset(){
    numParticles = 0;
    box.reset();
    energy = 0.0;
  }

  void grow(const Vector3D<Real> &v){
    box.grow(v);
  }

  /*
   * This method is called when performing a reduction over 
   * BoundingBox's. It subsumes the bounding box of the 'other'
   * and accumulates its energy in its own. If a PE has no
   * particles, its contributions are not counted.
   */
  void grow(const BoundingBox &other){
    if(other.numParticles == 0) return;
    if(numParticles == 0){
      *this = other;
    }
    else{
      box.grow(other.box);
      numParticles += other.numParticles;
      energy += other.energy;
    }
  }

  void expand(Real pad){
    box.greater_corner = box.greater_corner*pad+box.greater_corner;
    box.lesser_corner = box.lesser_corner-pad*box.lesser_corner;
  }

  void pup(PUP::er &p){
    p | box;
    p | numParticles;
    p | energy;
  }

};

#if 0
/* 
 * NodeDescriptor: used to keep track of the number of particles
 * on a particular PE under an active node during domain decomposition. 
 * It also tracks the largest and smallest keys of a node's particles
 * on a PE.
 */
struct NodeDescriptor {
  int numParticles;
  NodeDescriptor() {
    numParticles = 0;
  }

  NodeDescriptor(int np) : numParticles(np) {}

  /*
   * This method is called during the histogram reduction. Each PE
   * contributes a list of these descriptors, one for each active 
   * node that is being considered for partitioning. Through the reduction,
   * the total number of particles under an active node across all
   * PEs is obtained, as is the range of their keys. If a PE holds no
   * particles under a given node, its contribution is not considered.
   */
  void grow(const NodeDescriptor &other){
    numParticles += other.numParticles;
  }
};
#endif

class TreePiece;
struct ParticleMsg;
struct ForceData;
template<typename T> class Node;

/*
 * TreePieceDescriptor:
 * Used by the DataManager to keep track of the particles submitted by 
 * the TreePieces on its PE. It uses this information in the tree building
 * process to figure out which nodes are completely local to it (i.e. which
 * nodes have all of their particles on this PE). 
 */
struct TreePieceDescriptor {
  CkVec<ParticleMsg*> *vec;
  TreePiece *owner;
  int index;
  int numParticles;
  Node<ForceData> *root;

#ifndef SPLASH_COMPATIBLE
  int bucketStartIdx;
  int bucketEndIdx;
#endif

  TreePieceDescriptor() : 
    vec(NULL), owner(NULL), numParticles(0), index(-1)
  {
  }

  TreePieceDescriptor(CkVec<ParticleMsg*> *v, int np, TreePiece *o, int i) : 
    vec(v), owner(o), index(i), numParticles(np)
  {
  }

  TreePieceDescriptor(int i) : 
    vec(NULL), owner(NULL), numParticles(0), index(i)
  {
  }

  bool operator<=(const TreePieceDescriptor &t) const {
    return index <= t.index;
  }
  bool operator>=(const TreePieceDescriptor &t) const {
    return index >= t.index;
  }

  bool operator>=(const int &t) const {
    return index >= t;
  }
};

/*
 * ForceData:
 * Has fields to store the moments of a tree node. This includes the 
 * center of mass, total mass and radius of node.
 */
struct ForceData {
  OrientedBox<double> box;
  MultipoleMoments moments;

  ForceData() : box(), moments() 
  {
  }
};

/*
 * DtReductionStruct:
 * Used to tabulate basic statistics such as total number of interactions
 * and number of opening criterion calls.
 */
struct DtReductionStruct {
  CmiUInt8 pnInteractions;
  CmiUInt8 ppInteractions;
  CmiUInt8 openCrit;
#if 0
  bool haveNaN;
#endif

  DtReductionStruct &operator+=(const DtReductionStruct &other){
    pnInteractions += other.pnInteractions;
    ppInteractions += other.ppInteractions;
    openCrit += other.openCrit;
#if 0
    if(other.haveNaN) haveNaN = true;
#endif
    return *this;
  }
};

struct PointerContainer {
  void *ptr;
  PointerContainer() : ptr(NULL) {}
  PointerContainer(void *p) : ptr(p) {}
};

PUPbytes(PointerContainer);

#endif
