#ifndef __TREE_NODE_H__
#define __TREE_NODE_H__

/*
 * CharmBH: Node.h
 * Tree node data structure. Used in spatial partitioning of particles.
 * Each node has a key, type, depth, pointer to the particles underneath
 * it (valid only for "Internal" nodes, i.e. those nodes for which all
 * enclosed particles are hosted on the same PE as the node itself.), 
 * and child nodes. 
 * The Node data structure also tracks the owners of a node, i.e. 
 * which TreePieces have particles under that node. So, for example,
 * the root is owned collectively by all the TreePieces that have particles.
 */

#include "MultipoleMoments.h"
#include "OrientedBox.h"
#include "Particle.h"
#include "util.h"
#include "Descriptor.h"

#include <stack>
using namespace std;

// Type of node
enum NodeType {
  Invalid = 0,

  // All particles under the node are on this PE
  Internal,
  // This is a leaf
  Bucket,
  // This is a leaf with no particles in it
  EmptyBucket,

  /* 
   * Some of the particles underneath this node are
   * owned by TreePieces that are not on this PE
  */
  Boundary,

  /*
   * None of the particles underneath this node are owned 
   * by TreePieces on this PE, although they may have been
   * fetched from a remote source.
   */
  Remote,
  // Same as above, but this is a leaf
  RemoteBucket,
  // Same as above, but this is a leaf with no particles in it
  RemoteEmptyBucket
};

#if 0
extern string NodeTypeString[];
#endif

/*
 * These are "core" fields that would be present in any distributed memory 
 * tree data structure.
 * The part that can be customized for each algorithm is the 'data' field
 * of the Node data structure.
 */
struct NodeCore {
  Key key;
  NodeType type;

  int depth;
  // Pointer to particles underneath this node (valid only if Internal)
  Particle *particleStart;
  // Number of particles underneath this node
  int numParticles;
  // Number of children nodes 
  int numChildren;

  // TreePieces [ownerStart,ownerEnd] own this node
  int ownerStart;
  int ownerEnd;

  /* 
   * Was this node fetched from a remote source during traversal?
   * This check is needed to distinguish those nodes that were allocated
   * during the building of tree, from those that were received
   * as responses to requests for remote data. The distinction is important
   * because the former are allocated as arrays of size BRANCH_FACTOR (except
   * for the root, which is always allocated as a lone node), whereas the
   * latter are received as messages containing entire subtrees of nodes.
   */
  bool cached;

  NodeCore(){
    particleStart = NULL;
    numParticles = 0;
    numChildren = 0;
    type = Invalid;
    depth = -1;
    ownerStart = -1;
    ownerEnd = -1;
    cached = false;
  }

  NodeCore(Key k, int d, Particle *p, int np) : 
    key(k), depth(d), 
    particleStart(p), numParticles(np),
    type(Invalid), numChildren(0),
    ownerStart(-1), ownerEnd(-1), cached(false)
  {
  }
};

/*
 * Actual Node data structure. Has a core (described above)
 * and pointer to children (all BRANCH_FACTOR of which are allocated
 * in a single array) and the node's parent (NULL for root). 
 */
template<typename T>
class Node {

  NodeCore core;
  // all children are stored together
  Node<T> *children;
  Node<T> *parent;

  /* 
   * Used during tree building to ascertain whether the moments of all
   * of this node's children have been computed (or received). Then,
   * the moments of this node may be computed.
   */
  int numChildrenMomentsReady;

  public:
  /*
   * Data specific to the kind of tree we have built.
   * For a decomposition tree, this template would be instantiated 
   * with the NodeDescriptor type, whereas for the tree used during
   * traversal and force computation, it is of type ForceData, and
   * contains the MultipoleMoments of the node.
   */
  T data;

  Node(Key k, int d, Particle *p, int np, Node<T> *par=NULL) : 
    core(k,d,p,np), parent(par), children(NULL), data(),
    numChildrenMomentsReady(0)
  {
  }

  Node() : 
    core(), parent(NULL), children(NULL), data(),
    numChildrenMomentsReady(0)
  {
  }

  // What is the key of the first particle that can be enclosed by this node?
  static Key getParticleLevelKey(const Node<T> *node){
    Key k = node->getKey();
    int depth = node->getDepth();

    return (k<<(TREE_KEY_BITS-(LOG_BRANCH_FACTOR*depth+1)));
  }

  // ditto for last particle enclosed.
  static Key getLastParticleLevelKey(const Node<T> *node){
    Key k = node->getKey();
    int depth = node->getDepth();

    int nshift = TREE_KEY_BITS-(LOG_BRANCH_FACTOR*depth+1);
    k <<= nshift;
    Key mask = (Key(1) << nshift);
    mask--;
    k |= mask;
    return k;
  }

  // Find depth of node from position of prepended '1' bit in node
  static int getDepthFromKey(Key k){
    int nUsedBits = mssb64_pos(k);
    return (nUsedBits/LOG_BRANCH_FACTOR);
  }

  /*
   * Various get/set methods.
   */

  int getNumChildren() const {
    return core.numChildren;
  }

  int setChildren(Node<T> *ch, int n) {
    children = ch;
    core.numChildren = n;
  }

  Node<T> *getChildren() const { return children; }
  Node<T> *getChild(int i) const { return children+i; }
  int getNumParticles() const { return core.numParticles; }
  Particle *getParticles() const { return core.particleStart; }
  Key getKey() const { return core.key; }
  int getDepth() const { return core.depth; }
  int getOwnerStart() const { return core.ownerStart; }
  int getOwnerEnd() const { return core.ownerEnd; }
  NodeType getType() const { return core.type; }
  Node<T> *getParent() const { return parent; }
  void setKey(Key &k){ core.key = k; }
  void setDepth(int d){ core.depth = d; }
  void setType(NodeType t){ core.type = t; }
  void setParent(Node<T> *par){ parent = par; }

  void setOwners(int lb, int ub){
    core.ownerStart = lb;
    core.ownerEnd = ub;
  }

  void setOwnerStart(int lb){ core.ownerStart = lb; }
  void setOwnerEnd(int ub){ core.ownerEnd = ub; }

  void setParticles(Particle *p, int n){
    core.particleStart = p;
    core.numParticles = n;
  }

  void childMomentsReady(){ numChildrenMomentsReady++; }
  int getNumChildrenMomentsReady() const { return numChildrenMomentsReady; }
  bool allChildrenMomentsReady() const { return numChildrenMomentsReady == core.numChildren; }

#if 0
  void getOwnershipFromChildren(){
    // set node type from children types
    int numChildren = core.numChildren;
    bool allChildrenExternal = true;
    for(int i = 0; i < numChildren; i++){
      NodeType c_type = children[i].getType();
      if(c_type != Remote 
          && c_type != RemoteBucket 
          && c_type != RemoteEmptyBucket){
        allChildrenExternal = false;
      }
      if(i == 0) continue;
      int diff = children[i].getOwnerStart()
        - children[i-1].getOwnerEnd(); 
      CkAssert(diff >= 0 && diff <= 1);
    }

    // since we invoke this method only when updating
    // the local tree after moments exchange, the
    // nodes we encounter here can only be Remote/RemoteBucket/RemoteEmptyBucket
    // or Boundary
    if(allChildrenExternal) setType(Remote);
    else setType(Boundary);

    core.ownerStart = children[0].getOwnerStart();
    core.ownerEnd = children[numChildren-1].getOwnerEnd();

  }
#endif

  /*
   * refine() is invoked on a node that is to be partitioned. 
   * In this method, the node's children are allocated and 
   * initialized. Also, the particles of the parent node are
   * distributed among the children, as decided by the keys
   * of the children. 
   */
  void refine(){
    CkAssert(core.numChildren == 0);
    int numRankBits = LOG_BRANCH_FACTOR;
    int depth = getDepth();
    CkAssert(depth >= 0);
    CkAssert(depth < ((TREE_KEY_BITS-1)/numRankBits));

    // Allocate children
    core.numChildren = BRANCH_FACTOR;
    children = new Node<T>[BRANCH_FACTOR];

    Key myKey = getKey();

    int start = 0;
    int end = getNumParticles();
    Particle *particles = getParticles();
    int splitters[BRANCH_FACTOR+1];
    //splitters.resize(BRANCH_FACTOR+1);

    Key childKey = (myKey << numRankBits);
    int childDepth = depth+1;

    /*
     * findSplitters() distributes the parent's particles among the 
     * children. This is done by considering the Key of each child
     * and finding the first of the parent's particles whose Key is
     * geq to that of the child. This gives the starting point for the
     * particles of the child. The end is given by the start of the next
     * child. We are able to do this in-place partitioning because 
     * particles are assumed to have been sorted beforehand.
     */
    findSplitters(particles,start,end,splitters,childKey,childDepth);

    for(int i = 0; i < BRANCH_FACTOR; i++){
      initChild(i,splitters,childKey,childDepth);
      childKey++;
    }
  }

  // Initialize various fields of child based on those of parent.
  void initChild(int i, int *splitters, Key childKey, int childDepth){
    Node<T> *child = children+i;
    // The splitters array was filled in by findSplitters()
    int childPartStart = splitters[i]; 
    int childPartEnd = splitters[i+1]; 
    int childNumParticles = childPartEnd-childPartStart;

    // set child's key and depth
    child->setKey(childKey);
    child->setDepth(childDepth);
    child->setParent(this);

    // findSplitters distributed particles over
    // different children
    child->core.particleStart = core.particleStart+childPartStart;
    child->core.numParticles = childNumParticles;
    child->core.numChildren = 0;
  }

  /* 
   * serialize():
   * Called when the subtree under a node has to be communicated to a 
   * requestor DM. The method is recursive. Explanations of the code are
   * provided below.
   */
  void serialize(Node<T> *placeInBuf, Node<T> *&emptyBuf, int subtreeDepth){
    /* 
     * serialize() works by having each node copy its children into
     * the serialization buffer, and invoking serialize() on them
     * recursively. The procedure terminates when the entire subtree 
     * up to a certain depth under the original node has been placed
     * in the buffer, or if a node has no children.
     * 'emptybuf' points to that part of the array that hasn't yet 
     * been filled with nodes. This is where the children of this 
     * node will be placed. 'placeInBuf' is the pointer to the 
     * serialized copy of this node.
     */
    CmiUInt8 childOffset = emptyBuf-placeInBuf;
    // We have reached the depth cutoff for the subtree
    if(subtreeDepth == 0){
      // if subtreeDepth is 0 and placeInBuf is NULL, then
      // we have been asked to place all of the node's
      // successors up to depth 0, but not the root itself.
      // this is a nonsensical request.
      CkAssert(placeInBuf != NULL);
      // my children are not to be included
      placeInBuf->setChildren((Node<T>*)childOffset,0);
      return;
    }

    // PlaceInBuf is NULL for the root, since we don't want
    // to place the root itself in the serialization buffer, 
    // only the subtree under it, starting from its children.
    if(placeInBuf != NULL){
      // We exclude the root from the buffer 
      placeInBuf->setChildren((Node<T>*)childOffset,getNumChildren());
    }

    // Copy children into empty portion of buffer
    memcpy(emptyBuf,getChildren(),sizeof(Node<T>)*getNumChildren());
    // Set the parent pointers (in integer offsets) of your children
    // We don't do this for the serialized children of the subtree root, 
    // since they are being transported to a different PE.
    if(placeInBuf != NULL){
      Node<T> *childInBuf = emptyBuf;
      for(int i = 0; i < getNumChildren(); i++){
        childInBuf->setParent((Node<T>*)(placeInBuf-childInBuf));
        childInBuf++;
      }
    }
    // This is where I put the serialized copy of my first child
    placeInBuf = emptyBuf;
    // we have used up some of the empty space in the buffer
    emptyBuf += getNumChildren();

    // invoke serialize on each of the children in turn
    for(int i = 0; i < getNumChildren(); i++){
      getChild(i)->serialize(placeInBuf,emptyBuf,subtreeDepth-1);
      // move placeInBuf to point to the serialized copy of the next
      // child; emptyBuf is updated by serialize() recursively
      placeInBuf++;
    }
  }

  /*
   * This is the reverse of the serialize() procedure; it is invoked
   * on the DM (PE) that previously sent a request to the owner for 
   * (the subtree under) a node. Its purpose is to translate the 
   * integer offsets encoding parent-child relationships into normal
   * pointers so that these nodes can be traversed as a pointer-linked
   * tree data structure.
   * Here, 'start' points to the buffer of serialized nodes, and 'nn' 
   * is the number of nodes in that buffer. The method is invoked on
   * the parent on the requesting DM, whose children were requested 
   * in the first place. 
   */
  void deserialize(Node<T> *start, int nn){
    /*
     * Set the children of this node to be the first BRANCH_FACTOR 
       nodes in the buffer. This works because on the owner PE, the
       corresponding root of this subtree would have copied its 
       children into the serialization buffer first.
     */
    setChildren(start,BRANCH_FACTOR); 
    Node<T> *buf = start;
    // Go through all the nodes in the buffer and update their pointers
    for(int i = 0; i < nn; i++){
      // Set this node to be the parent of the first BRANCH_FACTOR nodes
      if(i < BRANCH_FACTOR){
        buf->setParent(this);
      }
      else{
        CmiUInt8 parentOffset = (CmiUInt8)(buf->getParent());
        // Add base address of this buffer to the relative offset
        // that was saved during serialization. This gives the absolute
        // address of the parent.
        buf->setParent(buf+parentOffset);
      }
      // Similar relative-to-absolute address conversion for pointers to
      // children 
      CmiUInt8 childrenOffset = (CmiUInt8)(buf->getChildren());
      buf->setChildren(buf+childrenOffset,buf->getNumChildren());

      // There are no particles under this remotely fetched node yet.
      buf->setParticles(NULL,0);
      // Translate type from what it was on the owner PE. e.g. a node that
      // was Internal on the owner PE becomes a Remote here, a Bucket becomes
      // a RemoteBucket, etc.
      buf->setType(makeRemote(buf->getType()));
      // This node has been fetched from a remote location, and is part of an
      // array that includes not just its siblings, other nodes too.
      buf->setCached();

      // Move to next node in serialized buffer.
      buf++;
    }
  }

  void setCached(){ core.cached = true; }
  bool isCached(){ return core.cached; }

  /*
   * Translate the type of a remotely fetched node.
   */
  static NodeType makeRemote(NodeType type){
    switch(type){
      case Boundary : 
      case Internal : return Remote;
      case Bucket : return RemoteBucket;
      case EmptyBucket : return RemoteEmptyBucket;

      case Remote :
      case RemoteBucket :
      case RemoteEmptyBucket : return type;

      default: CkAbort("bad type!\n");
    }
  }

  /*
   * Used during tree building. Calculate the moments of a node
   * from the moments of its children. 
   */
  void getMomentsFromChildren(){
    OrientedBox<Real> &bb = data.box;
    Real &mass = data.moments.totalMass;
    Vector3D<Real> &cm = data.moments.cm;

    for(Node<ForceData> *c = getChildren(); c != getChildren()+getNumChildren(); c++){
      MultipoleMoments &childMoments = c->data.moments;
      // My mass is the sum of the mass of each child
      mass += childMoments.totalMass;
      // My center of mass is the weighted average of the centers of mass 
      // of my children
      cm += childMoments.totalMass*childMoments.cm;
      // Grow my bounding box to enclose that of my child.
      bb.grow(c->data.box);
    }
    if(mass > 0.0) cm /= mass;

    // Radius: distance between center of mass and the corner that is
    // farthest from it.
    Vector3D<Real> delta1 = data.moments.cm - data.box.lesser_corner;	
    Vector3D<Real> delta2 = data.box.greater_corner - data.moments.cm;
    delta1.x = (delta1.x > delta2.x ? delta1.x : delta2.x);
    delta1.y = (delta1.y > delta2.y ? delta1.y : delta2.y);
    delta1.z = (delta1.z > delta2.z ? delta1.z : delta2.z);
    data.moments.rsq = delta1.lengthSquared();
  }

  /*
   * If this node is a leaf, obtain its moments from 
   * the particles that it encloses.
   */
  void getMomentsFromParticles(){
    OrientedBox<Real> &bb = data.box;
    Real &mass = data.moments.totalMass;
    Vector3D<Real> &cm = data.moments.cm;

    Particle *p;
    for(p = getParticles(); p != getParticles()+getNumParticles(); p++) {
      mass += p->mass;
      // Center of mass is weighted avg of particles' centers of mass
      cm += p->mass*p->position;
      bb.grow(p->position);
    }
    if(mass > 0.0) cm /= mass;

    Real d;
    // Radius: distance between center of mass and particle farthest from it
    data.moments.rsq = 0;
    for(p = getParticles(); p != getParticles()+getNumParticles(); p++) {
      d = (cm - p->position).lengthSquared();
      if(d > data.moments.rsq) data.moments.rsq = d;
    }

  }
};

/*
 * When a DM requests the information of a Remote node during tree building, 
 * it only needs certain fields, as embodied by this structure.
 */
struct MomentsExchangeStruct {
  MultipoleMoments moments;
  OrientedBox<Real> box;
  Key key;
  NodeType type;

  MomentsExchangeStruct() {}
  MomentsExchangeStruct(const Node<ForceData> &node):
    moments(node.data.moments), box(node.data.box),
    key(node.getKey()), type(node.getType()) {}
};
PUPbytes(MomentsExchangeStruct)


#endif
