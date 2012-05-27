struct TreeRootContainer {
  Node *root;
  int pe;

  TreeRootContainer(Node *r, int pe_) : 
    root(r),
    pe(pe_)
  {}
};

/*
 * BoundingBox:
 * Used to calculate the bounding box of all particles in the
 * simulation universe. It also keeps track of particle energy,
 * to see whether there is a drift in the total system energy 
 * through the simulation.
 */
class BoundingBox {
  OrientedBox<double> box;
  int numParticles;
  Real pe;
  Real ke;
  Real mass;

  public:
  
  BoundingBox(){
    reset();
  }

  void reset(){
    numParticles = 0;
    box.reset();
    pe = 0.0;
    ke = 0.0;
    mass = 0.0;
  }

  void grow(Vector3D<Real> v){
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
      pe += other.pe;
      ke += other.ke;
      mass += other.mass;
    }
  }

  void expand(Real pad){
    box.greater_corner = box.greater_corner*pad+box.greater_corner;
    box.lesser_corner = box.lesser_corner-pad*box.lesser_corner;
  }

  void pup(PUPer p){
    p|box;
    p|numParticles;
    p|pe;
    p|ke;
    p|mass;
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

  void copy(ForceData *other){
    box = other->box;
    moments = other->moments;
  }

  void reset(){
    *this = ForceData();
  }
};


/*
 * When a DM requests the information of a Remote node during tree building, 
 * it only needs certain fields, as embodied by this structure.
 */
class MomentsExchangeStruct {
  public MultipoleMoments moments;
  public OrientedBox<double> box;
  public Key key;
  public NodeType type;

  MomentsExchangeStruct() {}
  MomentsExchangeStruct(const Node &node){
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

class ParticleRequestTable {
  private:
  std::map<Key,ParticleRequest> table;

  public:
  void add(){
    // FIXME - flesh out
  }

  void free(){
    std::map<ParticleRequest>::iterator it;
    for(it = table.begin(); it != table.end(); it++){
      ParticleRequest &request = *it;
      CkAssert(request.sent);
      CkAssert(request.data != NULL);
      CkAssert(request.requestors.length() == 0);
      CkAssert(request.msg != NULL);

      CkAssert(request.parentCached == request.parent->isCached());

      if(!request.parentCached){
        delete request.parent;
      }

      delete request.msg;
    }
  }
  table.clear();
};

class NodeRequestTable {
  private:
  std::map<Key,NodeRequest> table;

  public:
  void free(){
    std::map<Key,NodeRequest>::iterator it;
    for(it = table.begin(); it != table.end(); it++){
      NodeRequest &request = *it;
      CkAssert(request.sent);
      CkAssert(request.data != NULL);
      CkAssert(request.requestors.length() == 0);
      CkAssert(request.msg != NULL);

      CkAssert(request.parentCached == request.parent->isCached());

      if(!request.parentCached){
        delete request.parent;
      }

      delete request.msg;
    }
    table.clear();
  }
};

template<typename T>
class HackArray {
  private:
  CkVec<T> array;

  HackArray() {} 
  HackArray(T *start, int n){
    array.resize(n);
    memcpy(&array[0],start,n*sizeof(T));
  }

  T *get(int i){
    return &array[i];
  }

  void add(T &elem){
    array.push_back(elem);
  }

  int &length(){
    return array.length();
  }

  void sort(){
    array.quickSort();
  }
};

typedef HackArray<Particle> ParticleArray;

struct PossiblySplitNode {
  Node *node;
  bool split;

  PossiblySplitNode() : 
    node(NULL),
    split(false)
  {}

  PossiblySplitNode(Node *n) : 
    node(n),
    split(false)
  {}
};

typedef HackArray<PossiblySplitNode> PossiblySplitNodeArray;

struct ActiveBinInfo{ 
  PossiblySplitNodeArray *oldvec;
  PossiblySplitNodeArray *newvec;
  CkVec<int> counts;

  ActiveBinInfo(){
    oldvec = new PossiblySplitNodeArray();
    newvec = new PossiblySplitNodeArray();
  }

  ~ActiveBinInfo(){
    delete oldvec;
    delete newvec;
  }

  /*
   * Add a new node to be considered for partitioning. 
   * Used to inititate the decomposition/tree building 
   * process with the root of the tree.
   */
  void addNewNode(Node *node){
    newvec->add(PossiblySplitNode(node));
    counts.push_back(node->getNumParticles());
  }

  /*
   * This method is invoked to apply the partitioning decisions
   * on the nodes in oldvec. 
   */
  void processRefine(CkVec<int> *binsToRefine){
    // For each index specfied in the list of refinements
    for(int i = 0; i < binsToRefine.size(); i++){
      // Get the current index of the node to refine
      int bin = (*binsToRefine)[i];
      PossiblySplitNode *possiblySplitNode = oldvec->get(bin);
      // Mark this node as refined.
      possiblySplitNode->split = true;
      refine(possiblySplitNode->node);
    }

  }

  virtual void refine(Node *node){
    node->refine(&counts,newvec,globalParams.decompLevels);
  }

  /*
   * These are used when sending the histogram contributions
   * of this PE to the partitioning reduction.
   */
  int getNumCounts(){
    return counts.length();
  }

  int getCount(int i){
    return counts[i];
  }

  /*
   * Swap oldvec and newvec, and clear counts;
   * also make space for new children by clearing
   * newvec.
   */
  void reset(){
    PossiblySplitNodeArray *tmp;
    tmp = oldvec;
    oldvec = newvec;
    newvec = tmp;
    newvec->length() = 0;

    counts.length() = 0;
  }

  PossiblySplitNodeArray *getActive(){
    return oldvec;
  }
};

class FlushParticlesStruct {
  public:
  CkVec<Key> *retractSites;
  int retractIndex;

  FlushParticlesStruct(CkVec<Key> *sites){
    retractIndex = 0;
    retractSites = sites;
  }

  bool doRetract(Key key){
    return (retractIndex < retractSites.length()
        && retractSites[retractIndex] == key);
  }

  void advance(){
    retractIndex++;
  }
};

/*
 * TreePieceDescriptor:
 * Used by the DataManager to keep track of the particles submitted by 
 * the TreePieces on its PE. It uses this information in the tree building
 * process to figure out which nodes are completely local to it (i.e. which
 * nodes have all of their particles on this PE). 
 */
class TreePieceDescriptor {
  CkVec<ParticleMsg*> *vec;
  TreePiece *owner;
  int index;
  int numParticles;
  Node *root;

  int bucketStartIdx;
  int bucketEndIdx;

  TreePieceDescriptor(){
    vec = NULL;
    owner = NULL;
    numParticles = 0; 
    index = -1;
  }

  TreePieceDescriptor(CkVec<ParticleMsg*> *v, int np, TreePiece o, int i){ 
    vec = v;
    owner = o;
    index = i;
    numParticles = np;
  }

  TreePieceDescriptor(int i){
    vec = NULL;
    owner = NULL;
    numParticles = 0;
    index = i;
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


typedef HackArray<TreePieceDescriptor> TreePieceDescriptorArray;
class TreePieceCounter : public CkLocIterator {            
  public:
  int count;
  int numParticles;
  TreePieceDescriptorArray submittedParticles;

  TreePieceCounter() { }                      

  void addLocation(CkLocation &loc) {
    const int *indexData = loc.getIndex().data();
    TreePiece *tp = treePieceProxy[indexData[0]].ckLocal();
    int np = tp->getNumParticles();
    submittedParticles.add(TreePieceDescriptor(tp->getBufferedParticleMsgs(), np, tp, indexData[0]));
    numParticles += np;
    count++;
  }

  void reset() {
    numParticles = 0;
    count = 0;
    submittedParticles.length() = 0;
  }

  void sortTreePieceDescriptors(){
    submittedParticles.sort();
  }

  TreePieceDescriptor *get(int i){
    return submittedParticles.get(i);
  }
};
