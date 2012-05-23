/*
 * TreePieceDescriptor:
 * Used by the DataManager to keep track of the particles submitted by 
 * the TreePieces on its PE. It uses this information in the tree building
 * process to figure out which nodes are completely local to it (i.e. which
 * nodes have all of their particles on this PE). 
 */
class TreePieceDescriptor {
  CkVec<ParticleMsg> vec;
  TreePiece owner;
  int index;
  int numParticles;
  Node<ForceData> root;

  int bucketStartIdx;
  int bucketEndIdx;

  TreePieceDescriptor(){
    vec = null;
    owner = null;
    numParticles = 0; 
    index = -1;
  }

  TreePieceDescriptor(CkVec<ParticleMsg> v, int np, TreePiece o, int i){ 
    vec = v;
    owner = o;
    index = i;
    numParticles = np;
  }

  TreePieceDescriptor(int i){
    vec = null;
    owner = null;
    numParticles = 0;
    index = i;
  }

  // FIXME - what to do with these?
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

