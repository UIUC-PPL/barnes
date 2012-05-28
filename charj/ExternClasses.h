struct TreeRootContainer {
  Node *root;
  int pe;

  TreeRootContainer(Node *r, int pe_) : 
    root(r),
    pe(pe_)
  {}
};

class ExternalParticle {
  Vector3D<Real> position;
  Real mass;

  void pup(PUP::er &p){
    p|position;
    p|mass;
  }
};

class Particle : public ExternalParticle {
  Key key;
  Vector3D<Real> velocity;
  Vector3D<Real> acceleration;
  Real potential;

  int order;

  boolean operator<=(const Particle &other) const { return key <= other.key; }
  boolean operator>=(const Particle &other) const { return key >= other.key; }
  boolean operator>=(const Key &k) const { return key >= k; }

  void pup(PUPer &p){
    ExternalParticle::pup(p);
    p|velocity;
    p|key;
    p|potential;
  }
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
  MomentsExchangeStruct(const Node *node){
    moments = node->data.moments;
    box = node->data.box;
    key = node->getKey()
    type = node->getType();
  }

  void pup(PUPer p){
    p|moments;
    p|box;
    p|key;
    p|type;
  }
};

struct NodeRequest {
  Key key;
  int replyTo;

  NodeRequest() {}
  NodeRequest(Key _key, int _replyTo){
    key = _key;
    replyTo = _replyTo;
  }
};
typedef HackArray<NodeRequest> NodeRequestArray;
const NodeRequest &make_NodeRequest(Key k, int reply);

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

  void pup(PUP::er &p){
    array.pup(p);
  }
};

typedef HackArray<Particle> ParticleArray;
typedef HackArray<ExternalParticle> ExternalParticleArray;

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

  void prepare(Node *root_, CkVec<Node*> *buckets){
    owner->prepare(root_,root,buckets,bucketStartIdx,bucketEndIdx);
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

  int &length(){
    return submittedParticles.length();
  }
};

struct NodeCacheLine {
  NodeReplyMsg *msg

  bool sent;
  bool haveData;
  CkVec<Requestor> requestors;

  Node *parent; 
  bool parentCached;

  NodeCacheLine() : 
    sent(false),
    haveData(false),
    msg(NULL),
    parent(NULL)
  {
  }

  void reset(){
    sent = false;
    haveData = false;
  }

  /*
   * When the subtree under a requested node is obtained on the PE, find the 
   * workers that requested it, and perform a traversal on this subtree
   * for each one. Note that the worker may be in the midst of a traversal for
   * a different bucket from the one that generated this request. Therefore,
   * we must save and restore its current bucket through get/setContext.
   */
  void deliver(){
    for(int i = 0; i < requestors.length(); i++){
      Requestor &req = requestors[i];
      CutoffWorker *worker = req.worker;

      // Save a pointer to the current bucket of the worker.
      Node *saveContext = worker->getContext();
      // Replace it with the one that generated the request
      // for this subtree
      worker->setContext(req.context);

      // Recall the kind of traversal the worker was performing
      // when it generated the request... 
      Traversal *traversal = req.traversal; 
      // And a pointer to its current book-keeping state
      State *state = req.state;
      // For each child of the parent whose children were requested...
      for(int j = 0; j < msg->getNumSubtreeRoots(); j++){
        // perform a top-down traversal of this subtree.
        traversal->topDownTraversal(msg->getSubtreeRoot(j),worker,state);
      }

      // Restore the current bucket of the worker
      worker->setContext(saveContext);
      // Check whether this traversal has completed.
      state->decrPending();
      if(state->complete()) worker->done();
    }
    requestors.clear();
  }


  void free(){
      CkAssert(sent);
      CkAssert(requestors.length() == 0);
      CkAssert(msg != NULL);

      CkAssert(parentCached == parent->isCached());

      // FIXME - might have to delete subtree 
      // in this cache line, depending on what
      // happens when msg is deleted
      if(!request.parentCached){
        delete request.parent;
      }

      delete request.msg;
  }
};

struct ParticleCacheLine {
  ParticleReplyMsg *msg

  bool sent;
  bool haveData;
  CkVec<Requestor> requestors;

  Node *parent; 
  bool parentCached;

  ParticleCacheLine() : 
    sent(false),
    haveData(false),
    msg(NULL),
    parent(NULL)
  {
  }

  void reset(){
    sent = false;
    haveData = false;
  }

  /*
   * When remote particles are received, we must supply them to all workers
   * on this PE that requested them. These workers will then use them for
   * computations that couldn't previously perform due to missing data.
   */
  void deliver(){
    // For each requestor of this set of particles
    for(int i = 0; i < requestors.length(); i++){
      Requestor &req = requestors[i];
      // Which worker requested these particles?
      CutoffWorker *worker = req.worker;

      /*
       * Workers maintain the current local bucket for which they
       * are performing the traversal (either remote or local). 
       * Currently, the worker could be in the middle of a traversal
       * for some other local bucket. We save the pointer to this 
       * local bucket here, and replace it with the one which had caused
       * the worker to request these particles in the first place.
       */
      Node *saveContext = worker->getContext();
      worker->setContext(req.context);

      // so that the worker may do some book-keeping, if required 
      worker->beforeParticleForces(parent->getKey());
      for(int j = 0; j < msg->getNumParticles(); j++){ 
        // supply each particle in the bucket to the worker
        worker->work(msg->getParticle(j));
      }
      // Tell the worker that it will receive no more particles from this
      // bucket
      worker->bucketDone(parent->getKey());

      // Restore the current bucket of the worker, in case
      // it is in the middle of some other traversal.
      worker->setContext(saveContext);

      // Update the corresponding state to reflect the fact that a previously
      // outstanding remote data request is now complete.
      State *state = req.state;
      state->decrPending();
      // If there are no more outstanding data requests for the 
      // traversal, we are done with it.
      if(state->complete()) worker->done();
    }

    // No more requestors of this bucket to keep track of 
    requestors.clear();
  }


  void free(){
      CkAssert(sent);
      CkAssert(requestors.length() == 0);
      CkAssert(msg != NULL);

      CkAssert(parentCached == parent->isCached());

      if(!parentCached){
        delete request.parent;
      }

      delete request.msg;
    }
  }
};

template <typename K, typename V>
class HackTable {
  std::map<K,V> table;

  V *get(K &k){
    std::map<K,V>::iterator it;
    it = table.find(k);
    if(it == table.end()){
      return NULL;
    }
    else{
      return &(it->second);
    }
  }

  V *put(K &k, const V &v = V()){
    std::pair<std::map<K,V>::iterator,bool> pr;
    pr = table.insert(make_pair(k,v));
    return &(pr.first->second);
  }

  void free(){
    std::map<K,V>::iterator it;
    for(it = table.begin(); it != table.end(); it++){
      V &value = *it;
      value.free();
    }
    table.clear();
  }
};

typedef HackTable<Key,NodeCacheLine> NodeCache;
typedef HackTable<Key,ParticleCacheLine> ParticleCache; 
typedef HackTable<Key,CkVec<int> > PendingMoments;

