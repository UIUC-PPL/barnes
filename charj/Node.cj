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

extern PossiblySplitNode;

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
class Node {
  Key key;
  NodeType type;

  int depth;
  // Pointer to particles underneath this node (valid only if Internal)
  Particle particles;
  int numParticles;

  Node parent;
  Node leftChild;
  Node rightChild;

  // TreePieces [ownerStart,ownerEnd) own this node
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
  boolean cached;


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
  ForceData data = new ForceData();

  Node(){
    type = Invalid;
    depth = -1;
    
    particles = null;
    numParticles = -1;

    ownerStart = -1;
    ownerEnd = -1;

    parent = null;
    leftChild = null;
    rightChild = null;

    cached = false;
    numChildrenMomentsReady = 0;
  }

  Node(Key k, int d, Particle p, int np, Node par) {
    key = k;
    depth = d;

    particles = p;
    numParticles = np;

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

  int getNumChildren() { return leftChild == null ? 0 : 2; }
  Node getLeftChild() { return leftChild; }
  Node getRightChild() { return rightChild; }
  void setLeftChild(Node node){ leftChild = node; }
  void setRightChild(Node node){ rightChild = node; }

  int getNumParticles() { return numParticles; }
  Particle getParticles() { return particles; }

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

  void setParticles(Particle p, int n){
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

    Key myKey = getKey();

    int start = 0;
    int end = start+getNumParticles();
    Particle p = getParticles();

    Array<int> splitters = new Array<int>([BRANCH_FACTOR+1]);

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

    leftChild = new Node();
    initChild(leftChild,this,splitters,childKey,childDepth);
    leftChild.refine(levels-1);
    childKey++;

    rightChild = new Node();
    initChild(rightChild,this,splitters,childKey,childDepth);
    rightChild.refine(levels-1);
    childKey++;

    delete splitters;
    
  }

  // this version is used during decomposition, 
  // in order to maintain a list of "active" nodes
  // for the next iteration of histogramming
  void refine(CkVec<int> counts, PossiblySplitNodeArray active){
    refine(counts,active,1);
  }

  void refine(CkVec<int> counts, PossiblySplitNodeArray active, int levels){
    if(levels == 0){
      active.add(this);
      counts.push_back(getNumParticles());
      return;
    }

    CkAssert(getNumChildren() == 0);
    int numRankBits = LOG_BRANCH_FACTOR;
    int depth = getDepth();
    CkAssert(depth >= 0);
    CkAssert(depth < ((TREE_KEY_BITS-1)/numRankBits));

    Key myKey = getKey();
    int start = 0;
    int end = start+getNumParticles();
    Particle p = getParticles();

    Array<int> splitters = new Array<int>([BRANCH_FACTOR+1]);

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

    leftChild = new Node();
    initChild(leftChild,this,splitters,childKey,childDepth);
    leftChild.refine(counts,active,levels-1);
    childKey++;

    rightChild = new Node();
    initChild(rightChild,this,splitters,childKey,childDepth);
    rightChild.refine(counts,active,levels-1);
    childKey++;

    delete splitters;
  }

  public NodeReplyMsg serialize(int subtreeDepth){
    NodeReplyMsg msg = new NodeReplyMsg(getKey());
    msg.roots[0] = getLeftChild().copy(subtreeDepth);
    msg.roots[1] = getRightChild().copy(subtreeDepth);
  }

  Node copy(int depth){
    Node root = new Node();
    root.copyFieldsFromNode(this);
    if(depth > 0 && getNumChildren() > 0){
      root.setLeftChild(getLeftChild().copy(depth-1));
      root.setRightChild(getRightChild().copy(depth-1));
    }
    return root;
  }

  void copyFieldsFromNode(Node other){
    key = other.getKey();
    type = other.getType();
    depth = other.getDepth();
    ownerStart = other.getOwnerStart();
    ownerEnd = other.getOwnerEnd();
    data.copy(other.getData());
  }

  void setCached(){ setCached(true); }
  void setCached(boolean val){ cached = val; }
  boolean isCached(){ return cached; }

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

  void deleteBeneath(){
    if(getNumChildren() > 0){
      getLeftChild().deleteRecursive();
      getRightChild().deleteRecursive();
    }
  }

  void deleteRecursive(){
    if(getNumChildren() > 0){
      getLeftChild().deleteRecursive();
      getRightChild().deleteRecursive();
    }
    delete this;
  }

  void reuseTree(){
    if(getType() == Invalid) return;
    reset();
    if(getNumChildren() == 0) return;
    for(int i = 0; i < getNumChildren(); i++){
      getChild(i)->reuseTree();
    }
  }

  void reset(){
    setType(Invalid);
    setParticles(null,-1);
    numChildrenMomentsReady = 0;
    data.reset();
  }

  boolean isLocal(){
    return (getType() == Internal ||
            getType() == Bucket   ||
            getType() == EmptyBucket);
            
  }

  boolean isBoundary(){
    return getType() == Boundary;
  }

  boolean isRemote(){
    return (getType() == Remote          ||
            getType() == RemoteBucket    ||
            getType() == RemoteEmptyBucket);
  }

  // FIXME - remember to set particles of node to null
  // before pup'ing
};


