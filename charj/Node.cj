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

/*
 * Actual Node data structure. Has a core (described above)
 * and pointer to children (all BRANCH_FACTOR of which are allocated
 * in a single array) and the node's parent (NULL for root). 
 */
class Node{
  Key key;
  NodeType type = Invalid;

  int depth = -1;
  // Pointer to particles underneath this node (valid only if Internal)
  Array<Particle> particles = null;
  int startParticleIdx;
  int numParticles;

  Array<Node> children = null; 
  int numChildren = 0;

  // TreePieces [ownerStart,ownerEnd) own this node
  int ownerStart = -1;
  int ownerEnd = -1;

  /* 
   * Was this node fetched from a remote source during traversal?
   * This check is needed to distinguish those nodes that were allocated
   * during the building of tree, from those that were received
   * as responses to requests for remote data. The distinction is important
   * because the former are allocated as arrays of size BRANCH_FACTOR (except
   * for the root, which is always allocated as a lone node), whereas the
   * latter are received as messages containing entire subtrees of nodes.
   */
  bool cached = false;

  Node parent;

  /* 
   * Used during tree building to ascertain whether the moments of all
   * of this node's children have been computed (or received). Then,
   * the moments of this node may be computed.
   */
  int numChildrenMomentsReady = 0;

  public:
  /*
   * Data specific to the kind of tree we have built.
   * For a decomposition tree, this template would be instantiated 
   * with the NodeDescriptor type, whereas for the tree used during
   * traversal and force computation, it is of type ForceData, and
   * contains the MultipoleMoments of the node.
   */
  ForceData data;

  Node(Key k, int d, Array<Particle> p, int start, int np, Node par) {
    key = k;
    depth = d;
    particles = p;
    numParticles = np;
    startParticleIdx = start;
    type = Invalid;
    parent = par;

    numChildrenMomentsReady = 0;
  }

  // What is the key of the first particle that can be enclosed by this node?
  static Key getParticleLevelKey(Node node){
    Key k = node.getKey();
    int depth = node.getDepth();

    return (k<<(TREE_KEY_BITS-(LOG_BRANCH_FACTOR*depth+1)));
  }

  // ditto for last particle enclosed.
  static Key getLastParticleLevelKey(Node node){
    Key k = node.getKey();
    int depth = node.getDepth();

    int nshift = TREE_KEY_BITS-(LOG_BRANCH_FACTOR*depth+1);
    k <<= nshift;
    Key mask = (((Key)1) << nshift);
    mask--;
    k |= mask;
    return k;
  }

  // Find depth of node from position of prepended '1' bit in node
  static int getDepthFromKey(Key k){
    int nUsedBits = mssb64_pos(k);
    return (nUsedBits/LOG_BRANCH_FACTOR);
  }

  static int completeTreeSize(int levels){
    return numLeaves(levels+1)-1;
  }

  static int numLeaves(int levels){
    return (1<<(LOG_BRANCH_FACTOR*levels));
  }
  /*
   * Various get/set methods.
   */

  int getNumChildren() {
    return numChildren;
  }

  void setChildren(Array<Node> ch, int n) {
    children = ch;
    numChildren = n;
  }

  Array<Node> getChildren() { return children; }
  Node getChild(int i) { return children[i]; }
  int getNumParticles() { return numParticles; }
  Array<Particle> getParticle(int i) { return particles[i]; }
  Key getKey() { return key; }
  int getDepth() { return depth; }
  int getOwnerStart() { return ownerStart; }
  int getOwnerEnd() { return ownerEnd; }
  NodeType getType() { return type; }
  bool isInternal() {
    return (type == EmptyBucket || type == Bucket || type == Internal);
  }
  Node getParent() { return parent; }
  void setKey(Key k){ key = k; }
  void setDepth(int d){ depth = d; }
  void setType(NodeType t){ type = t; }
  void setParent(Node par){ parent = par; }

  void setOwners(int lb, int ub){
    ownerStart = lb;
    ownerEnd = ub;
  }

  void setOwnerStart(int lb){ ownerStart = lb; }
  void setOwnerEnd(int ub){ ownerEnd = ub; }

  void setParticles(Array<Particle> p, int start, int n){
    startParticleIdx = start
    particles = p;
    numParticles = n;
  }

  void setNumParticles(int n){
    numParticles = n;
  }

  void childMomentsReady(){ numChildrenMomentsReady++; }
  int getNumChildrenMomentsReady() { return numChildrenMomentsReady; }
  bool allChildrenMomentsReady() { return numChildrenMomentsReady == getNumChildren(); }
  void setChildrenMomentsReady() { numChildrenMomentsReady = getNumChildren(); }

  /*
   * refine() is invoked on a node that is to be partitioned. 
   * In this method, the node's children are allocated and 
   * initialized. Also, the particles of the parent node are
   * distributed among the children, as decided by the keys
   * of the children. 
   */
    
  void refine(){
    refine(1);
  }

  void refine(int levels){
    if(levels == 0) return;

    CkAssert(getNumChildren() == 0);
    int numRankBits = LOG_BRANCH_FACTOR;
    int depth = getDepth();
    CkAssert(depth >= 0);
    CkAssert(depth < ((TREE_KEY_BITS-1)/numRankBits));

    // Allocate children
    numChildren = BRANCH_FACTOR;
    children = new Array<Node>(BRANCH_FACTOR);

    Key myKey = getKey();

    int start = 0;
    int end = getNumParticles();
    Array<Particle> p = getParticles();
    Array<int> splitters = new Array<int>(BRANCH_FACTOR+1);

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
    findSplitters(p,start,end,splitters,childKey,childDepth);

    for(int i = 0; i < BRANCH_FACTOR; i++){
      initChild(i,splitters,childKey,childDepth);
      getChild(i)->refine(levels-1);
      childKey++;
    }

    delete splitters;
    
  }

  // this version is used during decomposition, 
  // in order to maintain a list of "active" nodes
  // for the next iteration of histogramming
  void refine(CkVec<int> counts, CkVec<PossiblySplitNode > active){
    refine(counts,active,1);
  }

  void refine(CkVec<int> counts, CkVec<PossiblySplitNode > active, int levels){
    if(levels == 0){
      active.push_back(new PossiblySplitNode(this));
      counts.push_back(getNumParticles());
      return;
    }

    CkAssert(getNumChildren() == 0);
    int numRankBits = LOG_BRANCH_FACTOR;
    int depth = getDepth();
    CkAssert(depth >= 0);
    CkAssert(depth < ((TREE_KEY_BITS-1)/numRankBits));

    // Allocate children
    numChildren = BRANCH_FACTOR;
    children = new Array<Node>(BRANCH_FACTOR);

    Key myKey = getKey();

    int start = 0;
    int end = getNumParticles();
    Array<Particle> p = getParticles();
    Array<int> splitters = new Array<int>(BRANCH_FACTOR+1);

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
    findSplitters(p,start,end,splitters,childKey,childDepth);

    for(int i = 0; i < BRANCH_FACTOR; i++){
      initChild(i,splitters,childKey,childDepth);
      getChild(i).refine(counts,active,levels-1);
      childKey++;
    }

    delete splitters;
  }

  void reuseRefine(){
    int start = 0;
    int end = getNumParticles();
    Array<Particle> p = getParticles();

    int numRankBits = LOG_BRANCH_FACTOR;
    Key childKey = (getKey() << numRankBits);
    int childDepth = getDepth()+1;

    Array<int> splitters = new Array<int>(BRANCH_FACTOR+1);
    findSplitters(p,start,end,splitters,childKey,childDepth);
    Array<Node> child = getChildren();
    for(int i = 0; i < BRANCH_FACTOR; i++){
      start = splitters[i];
      end = splitters[i+1];
      child.setParticles(p,start,end-start);
      child++; 
    }
  }

  // Initialize various fields of child based on those of parent.
  void initChild(int i, Array<int> splitters, Key childKey, int childDepth){
    Node child = getChild(i);
    // The splitters array was filled in by findSplitters()
    int childPartStart = splitters[i]; 
    int childPartEnd = splitters[i+1]; 
    int childNumParticles = childPartEnd-childPartStart;

    // set child's key and depth
    child.setKey(childKey);
    child.setDepth(childDepth);
    child.setParent(this);

    // findSplitters distributed particles over
    // different children
    child.setParticles(getParticles(),childPartStart,childNumParticles);
    child.setChildren(null,0);
  }

  NodeReplyMsg serialize(int subtreeDepth){
    NodeReplyMsg msg = new NodeReplyMsg(getKey());
    for(int i = 0; i < BRANCH_FACTOR; i++){
      msg.roots[i] = getChild(i).copy(subtreeDepth);
    }
  }

  Node copy(int depth){
    Node root = new Node(this);
    if(depth > 0 && getNumChildren() > 0){
      root.setChildren(new Array<Node>(BRANCH_FACTOR),BRANCH_FACTOR);
      for(int i = 0; i < getNumChildren(); i++){
        root.getChild(i) = copy(getChild(i),depth-1);
      }
    }
    else{
      root.setChildren(null,0);
    }
    return root;
  }

  void setCached(){ setCached(true); }
  void setCached(boolean val){ cached = val; }
  bool isCached(){ return cached; }

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
    return Invalid;
  }

  /*
   * Used during tree building. Calculate the moments of a node
   * from the moments of its children. 
   */
  void getMomentsFromChildren(){
    OrientedBox<double> bb = data.box;
    Real mass = data.moments.totalMass;
    Vector3D<Real> &cm = data.moments.cm;

    for(int i = 0; i < getNumChildren(); i++){
      Node child = getChild(i);
      MultipoleMoments childMoments = child.data.moments;
      // My mass is the sum of the mass of each child
      mass += childMoments.totalMass;
      // My center of mass is the weighted average of the centers of mass 
      // of my children
      // FIXME += allowed?
      cm += childMoments.totalMass*childMoments.cm;
      // Grow my bounding box to enclose that of my child.
      bb.grow(child.data.box);
    }
    if(mass > 0.0) cm /= mass;

    // Radius: distance between center of mass and the corner that is
    // farthest from it.
    // FIXME 
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
    OrientedBox<double> bb = data.box;
    Real mass = data.moments.totalMass;
    Vector3D<Real> cm = data.moments.cm;

    Array<Particle> p = getParticles();
    for(int i = 0; i < getNumParticles(); i++){
      mass += p.mass;
      // Center of mass is weighted avg of particles' centers of mass
      cm += p.mass*p.position;
      bb.grow(p.position);
    }
    if(mass > 0.0) cm /= mass;

    Real d;
    // Radius: distance between center of mass and particle farthest from it
    data.moments.rsq = 0;
    for(p = getParticles(); p != getParticles()+getNumParticles(); p++) {
      d = (cm - position).lengthSquared();
      if(d > data.moments.rsq) data.moments.rsq = d;
    }

  }

  void deleteBeneath(){
    if(getNumChildren() == 0) return;
    for(int i = 0; i < getNumChildren(); i++){
      getChild(i).deleteBeneath();
    }
    delete getChildren();
    setChildren(null,-1,0);
  }

  void reuseTree(){
    if(getType() == Invalid) return;
    reset();
    if(getNumChildren() == 0) return;
    for(int i = 0; i < getNumChildren(); i++){
      getChild(i).reuseTree();
    }
  }

  void reset(){
    setType(Invalid);
    setParticles(null,-1,0);
    numChildrenMomentsReady = 0;
    data = T();
  }

};

/*
 * When a DM requests the information of a Remote node during tree building, 
 * it only needs certain fields, as embodied by this structure.
 */
class MomentsExchangeStruct {
  public MultipoleMoments moments = new MultipoleMoments();
  public OrientedBox<double> box = new OrientedBox<double>();
  public Key key;
  public NodeType type = Invalid;

  MomentsExchangeStruct(Node<ForceData> node){
    moments = node.data.moments;
    box = node.data.box;
    key = node.getKey()
    type = node.getType();
  }

  void pup(PUPer p){
    p|moments;
    p|box;
    p|key;
    p|type;
  }
};

class NodeRequest {
  Key key;
  int replyTo;

  NodeRequest() {}
  NodeRequest(Key _key, int _replyTo){
    key = _key;
    replyTo = _replyTo;
  }

  Key getKey() { return key; }
  int getReplyTo() { return replyTo; }
};
