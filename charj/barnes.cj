/*
  Proxies for tree pieces, data managers and main chare.
*/
readonly TreePiece@ treePieceProxy;
readonly DataManager@ dataManagerProxy;
readonly TreeMerger@ treeMergerProxy;
readonly Main@ mainProxy;
readonly Parameters globalParams;

// used while decomposing particles onto tree pieces
message ParticleMsg {
  ParticleArray part;
};

// To respond to remote request for particles 
message ParticleReplyMsg {
  Key key;
  ExternalParticleArray particles;

  ParticleReplyMsg(Key k, ExternalParticleArray p){
    key = k;
    particles = p;
  }
};

// To respond to remote request for subtree under node
message NodeReplyMsg {
  // key of parent of the subtrees sent 
  // via this message
  public Key key;
  public Array<Node> roots = new Array<Node>(BRANCH_FACTOR);

  public NodeReplyMsg(Key key_){
    key = key_;
  }
};

// To continuously reschedule local and remote work
message RescheduleMsg;

extern "C++" void setParameters(CkArgMsg *msg, Parameters &globalParams);
extern Parameters;

// Main chare initiates computation and ends it
mainchare Main {
  entry Main(CkArgMsg msg){
    /* Read input parameters from command line. Fill in defaults for
       unspecified parameters.
     */
    setParameters(msg);

    /* Create the DataManager group; there is one representative/member
       of the group on each PE.
     */
    treeMergerProxy = new TreeMerger@();
    dataManagerProxy = new DataManager@();

    CkArrayOptions opts(globalParams.numTreePieces);
    /* 
       We would like the tree piece array elements to be mapped to PEs
       in a round-robin fashion to begin with.
     */
    CProxy_RRMap myMap = CProxy_RRMap::ckNew();
    opts.setMap(myMap);
    /*
       Create the tree piece chare array.
     */
    treePieceProxy = new TreePiece@(opts);

    /*
       Set mainProxy readonly so that it is available
       on all other PEs after this method finishes.
     */
    mainProxy = thisProxy;

    /*
       Begin actual computation
     */
    thisProxy@commence();

    delete msg;
  }

  entry void commence(){
    double loadTime = CmiWallTimer();
    CkPrintf("[Main] load particles ... ");
    /* 
       Tell each PE to read the particles 
       from its portion of the input file
     */
    dataManagerProxy@loadTipsy(CkCallback(CkReductionTarget(Main,doneLoadParticles),thisProxy));
  }

  entry void doneLoadParticles(Array<Real> bb){
    //dataManagerProxy.loadParticles(CkCallbackResumeThread((void *&)redMsg));

    CkPrintf(" took %f s\n", CmiWallTimer()-loadTime);
    CkPrintf("bb l %f %f %f g %f %f %f\n", bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);

    //CkExit();
    //return;

    /*
       Each PE contributes the bounding box of the particles
       it read to a reduction; as a result, the bounding box
       of all particles in the simulated "universe" is obtained.
     */
    BoundingBox universe;
    //Real *data = (Real *) redMsg->getData();
    for(int i = 0; i < 3; i++){
      universe.box.lesser_corner[i] = bb[i];
      universe.box.greater_corner[i] = bb[i+3];
    }
    universe.pe = data[6];
    universe.ke = data[7];

    Real pad = 0.00001;
    /*
       Provide a little padding so that the particles closest to
       the extremities of the universe get correct keys.
     */
    universe.expand(pad);
    /*
       Tell all PEs to begin Oct decomposition of read particles.
     */

    CkStartQD(CkCallback(CkReductionTarget(DataManager,doneParticleFlush),0,dataManagerProxy));
    dataManagerProxy@decompose(universe);

  }

  void niceExit(){
    CkPrintf("[Main] graceful exit\n");
    CkExit();
  }

  void quiescence(){
    CkPrintf("[Main] Quiescence detected!\n");
    treePieceProxy.quiescence();
    dataManagerProxy.quiescence();
  }

  void quiescenceExit(){
    numQuiescenceRecvd++;
    if(numQuiescenceRecvd == 2){
      CkExit();
    }
  }


  Parameters params;
  int numQuiescenceRecvd;
};


  /*
    There is one DataManager per PE. It has several functions: 
    1. load a subset of the particles from an input file
    2. Oct-decompose particles onto TreePieces using a histogram-sort like algorithm
    3. Construct a PE-level tree from all the particles on the PE after decomposition
    4. Initiate traversals for all TreePieces on the PE after tree building
    5. Maintain a per-PE cache of data fetched from remote sources for TreePieces on the PE
    6. Update particle positions and velocities once traversals have completed.
  */
group DataManager {

  // data members
  int numRankBits;
  int numPesPerNode;
  double prevIterationStart;

  double phaseTime;

  ParticleArray myParticles;
  int myNumParticles;
  BoundingBox myBox;
  BoundingBox univBox;
  Real uside;

  int numTreePieces;
  int numLocalUsefulTreePieces;
  TreePieceCounter localTreePieces;

  int numMomentsRequested;
  int numMomentsReceived;

  int iteration;
  int decompIterations;
  ActiveBinInfo activeBins;
  // Root of particles on this PE
  // Used in both decomposition and traversal
  Node root;

  CkVec<Node> myBuckets;

  // I am done constructing the tree 
  // from particles present on this PE
  // However, in general, there will still be nodes whose
  // moments haven't yet been computed, since they might require 
  // remote information
  bool doneTreeBuild;
  std.map<Key,Node> nodeTable;

  PendingMoments pendingMoments;

  // I have processed the moment 
  // contributions from all other PEs, so that
  // the tree on this PE is now ready for 
  // traversal
  bool treeMomentsReady;
  // requests received from remote tree pieces
  // for node moments during tree building phase
  NodeRequestArray bufferedNodeRequests;
  NodeRequestArray bufferedParticleRequests;

  // requests made by local tree pieces for 
  // remote data during tree traversal
  NodeCache nodeCache;
  ParticleCache particleCache;

  // needed for skipping decomposition
  KeyArray treePieceKeys;
  boolean doSkipDecomposition;

  TreeSizeWorker treeSizeWorker = new TreeSizeWorker();
  Traversal fillTrav = new Traversal();
  int numTreePiecesDoneTraversals;
  int numTreePiecesDoneRemoteRequests;
  CacheStats nodeReqs;
  CacheStats partReqs;

  /* used only on PE 0 to track drift in energy */ 
  Real compareEnergy;

  // FIXME is this 64-bit?
  Array<long> numInteractions = new Array<long>([3]);

  entry DataManager(){
    numRankBits = LOG_BRANCH_FACTOR;
    iteration = 0;
    prevIterationStart = 0.0;
    root = null; 
    numPesPerNode = CkNumPes()/CkNumNodes();
    numLocalUsefulTreePieces = 0;
    doSkipDecomposition = false;
    init();
  }

  entry void loadTipsy(CkCallback cb){
    string filename(globalParams.filename);
    Tipsy.TipsyReader r = new Tipsy.TipsyReader(filename);
    if(!r.status()) {
      CkPrintf("[%d] : Fatal: Couldn't open tipsy file (%s)\n", CkMyPe(), filename);
      CkExit();
      return;
    }

    
    int nTotalParticles = getTipsyNBodies(r);
    int excess;
    int startParticle;

    int numLoaders = CkNumPes();

    myNumParticles = nTotalParticles / numLoaders;
    excess = nTotalParticles % numLoaders;
    startParticle = myNumParticles * CkMyPe();
    if(CkMyPe() < (int) excess) {
      myNumParticles++;
      startParticle += CkMyPe();
    }
    else {
      startParticle += excess;
    }

    // allocate an array for myParticles
    int nStore = myNumParticles;
    myParticles.resize(nStore);

    if(!r.seekParticleNum(startParticle)) {
      CkAbort("Couldn't seek to my particles!");
      return;
    }

    myBox.reset();
    myBox.pe = 0.0;
    myBox.ke = 0.0;

    for(int i = 0; i < myNumParticles; ++i) {
      Particle p = myParticles.get(i);
      setTipsyDarkParticle(p,r);
      p.order = i+startParticle; 
      p.potential = 0.0;

      myBox.grow(p.position);
      myBox.mass += p.mass;
      myBox.ke += p.mass*p.velocity.lengthSquared();
      myBox.pe = 0.0;
    }
    myBox.ke /= 2.0;
    contributeBoundingBox(cb);

    delete r;
  }

  /*
     Use the keys of the particles (as calculated above) 
     to decompose them onto tree pieces. Each tree piece
     gets a subvolume of the entire universe. This subvolume
     is represented in the tree by a particular node (internal
     or leaf). All particles within this subvolume fall beneath
     the node, so that the node's key is the common prefix of 
     all particles enclosed by the subvolume.

     Through the decomposition, we obtain several subsets of 
     particles, each of which constitutes a subvolume of the 
     universe. This subvolume corresponds to a tree piece, and
     the number of particles within it is guaranteed to be less
     than a user-specified threshold (ppc/particles per chare).

     The decomposition procedure is iterative, similar to
     histogram sort. In each iteration, a master PE (0) maintains
     a list of "active" nodes (leaves) which are to be partitioned.
     The keys of these nodes are broadcast to the workers (PE 0 is 
     both master and worker). The workers split the corresponding
     nodes into children and partition the particles previously in
     the active nodes to the appropriate children. They then contribute
     the number of particles they hold under each child to a reduction.
     The master uses the result of this reduction to find which nodes
     are to be split further: These are made active leaves for the 
     next iteration; this continues until there are no remaining active
     leaves.

     Decomposition begins at the root of the global tree, which 
     represents the bounding box enclosing all particles in the 
     simulated universe.
   */
  entry void decompose(BoundingBox universe){
    univBox = universe;

    uside = (univBox.box.greater_corner-univBox.box.lesser_corner).lengthSquared();
    uside *= 1.00002;

    if(CkMyPe() == 0){
      phaseTime = CmiWallTimer();
    }

    // Obtain key for each particle given the extents of the universe.
    hashParticleCoordinates(univBox.box);
    // Sort particles so that decomposition can be done in-place.
    myParticles.quickSort();

    if(CkMyPe()==0){
      Real deltaE = 0.0;
      // Check whether total energy remains (about) constant
      if(iteration == 1){
        // save this value so that we can compare
        // against it in future iterations. 
        compareEnergy = univBox.ke+univBox.pe;
      }
      else if(iteration > 1){
        // FIXME - disabled until we use custom reducer for bounding box
        /*
        deltaE = compareEnergy-(univBox.ke+univBox.pe);
        deltaE /= compareEnergy;
        if(deltaE < 0.0) deltaE = -deltaE;
        */
      }

      // Print statistics
      float memMB = (1.0*CmiMemoryUsage())/(1<<20);
      /*
      ostringstream oss; 
      CkPrintf("(%d) iteration %d rsize %f energy %f delE/E %f\n", 
          CkMyPe(),
          iteration,
          uside,
          univBox.ke+univBox.pe, deltaE);
      */

      CkPrintf("(%d) iteration %d rsize %f mem %.2f MB prev time %g s\n\n", CkMyPe(), iteration, uside, memMB, CmiWallTimer()-prevIterationStart);

      if(iteration > 1){
        // The energy should grow in magnitude
        // by less than a tenth of one per cent.
        //CkAssert(deltaE < 0.01);
      }

      prevIterationStart = CkWallTimer();
    }

    if(doSkipDecomposition){
      skipFlushParticles();
    }
    else{
      // How many particles do I hold?
      initHistogramParticles();
      // Send this count to the master PE
      sendHistogram();
    }
  }

  // executed on PE 0
  entry void receiveHistogram(Array<int> counts){
    int numRecvdBins = counts.size(0);
    CkVec<int> binsToRefine = new CkVec<int>;

    /* 
       Reserve enough space so that all active leaves 
       can be partitioned, if necessary
     */
    binsToRefine.reserve(numRecvdBins);
    binsToRefine.length() = 0;

    int particlesHistogrammed = 0;
    /*
       The newly created children of the active leaves
       are now the active leaves themselves. 
     */
    PossiblySplitNodeArray active = activeBins.getActive();
    CkAssert(numRecvdBins == active.length());

    Real thresh = (DECOMP_TOLERANCE*(Real)(globalParams.ppc));
    // Check which new active leaves need to be partitioned
    for(int i = 0; i < numRecvdBins; i++){
      PossiblySplitNode possiblySplitNode = active.get(i);
      Node node = possiblySplitNode.node;
      if((Real)counts[i] > thresh){
        // Need to refine this leaf (partition)
        binsToRefine.push_back(i);
        // By refining this node, we will remove one tree piece
        // and add BRANCH_FACTOR in its place.
        //CkPrintf("[%d] Split bin %d num %d thresh %.1f\n", CkMyPe(), i, counts[i], thresh);
        // IMPORTANT: we use ownerstart as a placeholder for 
        // global num. particles underneath 'node'
        node.setOwnerStart(-1);
      }
      else{
        // this node will not be refined; save the
        // global number of particles underneath it
        // we will use this count later, to decide
        // which nodes can be retracted.
        node.setOwnerStart(counts[i]);
      }

      particlesHistogrammed += counts[i];
    }

    // Do any active leaves need to be partitioned?
    if(binsToRefine.size() > 0) {
      // Yes, tell workers which ones  
      thisProxy@receiveSplitters(binsToRefine);
      decompIterations++;
    }
    else{
      // No more leaves to refine; send out particles to tree pieces
      double t = CmiWallTimer();
      CkPrintf("[0] Decomposition took %f s iterations %d\n", t-phaseTime, decompIterations);
      phaseTime = t;
      decompIterations = 0;
      doneDecomposition();
    }

    delete binsToRefine;
    // FIXME - do i have to delete counts?
    delete counts;
  }

  /*
     The master communicates the active leaves to split
     to the workers through this method. The activeBins
     data structure processes the list of leaves to refine,
     creates new children for them and obtains counts
     of the number of particles under each child. These
     children are then made the active leaves for the new
     iteration and the counts are contributed to a reduction
     in sendHistogram.
   */
  entry void receiveSplitters(CkVec<int> splitBins){

    // Process bins to refine
    activeBins.processRefine(splitBins);

    // Send counts of particles in active leaves to master PE
    sendHistogram();
  }


  /*
     This method is invoked on workers by the master, telling
     them that it is now OK to send the particles they are
     holding to the tree pieces that they are meant for.

     It also communicates the ranges of particles held  by each
     tree piece to every PE. 
   */
  entry void sendParticles(CkVec<Key> retractSites){
    // Flush particles to their owner tree pieces
    // this will also set the number of tree pieces used
    // after retracting
    FlushParticlesStruct fps = new FlushParticlesStruct(retractSites);
    numTreePieces = flushParticles(fps);
    CkAssert(fps.retractIndex == fps.retractSites.length());

    // Fail if not enough tree pieces are provided
    if(CkMyPe() == 0 && numTreePieces > globalParams.numTreePieces){
      CkPrintf("have %d treepieces need %d\n",globalParams.numTreePieces,numTreePieces);
      CkAbort("Need more tree pieces!\n");
    }
    else if(CkMyPe() == 0){
      CkPrintf("Used %d tree pieces\n", numTreePieces);
    }
  }

  entry void receiveMoments(MomentsExchangeStruct moments);

  entry void processSubmittedParticles();
  entry void doneStreaming();
  entry void startTraversal();

  entry void recvParticles(ParticleReplyMsg msg);
  entry void recvNode(NodeReplyMsg msg);
  entry void advance(CkReductionMsg msg);
  entry void traversalsDone(CmiUInt8 pnInter, CmiUInt8 ppInter, CmiUInt8 openCrit);

  entry void recvUnivBoundingBox(CkReductionMsg msg);

  entry void quiescence();
  entry void resumeFromLB();

  entry void doneRemoteRequests();

  entry void doneParticleFlush(){
    double t = CmiWallTimer();
    CkPrintf("[0] Flush took %f s\n", t-phaseTime);
    phaseTime = t;

    thisProxy@resumeProcessSubmittedParticles();
  }

  entry void resumeProcessSubmittedParticles(){
    // get the tree pieces (and their particles) on this PE
    senseTreePieces();

    // sort the local tree pieces by index (and hence lower/upper bound of particle range)
    localTreePieces.sortTreePieceDescriptors();
    myNumParticles = localTreePieces.numParticles;

    myParticles.length() = 0;
    myParticles.resize(myNumParticles);

    int offset = 0;
    numLocalUsefulTreePieces = 0;
    for(int i = 0; i < localTreePieces.count; i++){
      TreePieceDescriptor descr = localTreePieces.get(i);
      if(descr.index < numTreePieces) numLocalUsefulTreePieces++;

      offset += copyParticlesFromTreePieceDescriptor(myParticles.get(offset),descr);
    }

    myParticles.sort();

    // The first local TreePiece will always have buckets numbered from 0
    if(localTreePieces.count > 0) localTreePieces.get(0).bucketStartIdx = 0;

    buildTree(root,0,myNumParticles,0,numLocalUsefulTreePieces);
    doneTreeBuild = true;
  }

  entry void doneTreeBuilds();
  entry void resumeDoneNodeLevelMerge();
  entry void doneForces();
  entry void resumeFinishIteration();

  // sequential methods
  void copyMomentsToNode(Node<ForceData> *node, MomentsExchangeStruct mes){
    CkAssert(node.getKey() == mes.key);

    node.data.moments = mes.moments;
    node.data.box = mes.box;
    NodeType type = mes.type;
    node.setType(Node<ForceData>.makeRemote(type));
  }

  /*
     Obtain a 64-bit key for each particle. This gives us a cheap way to
     arrange particles within the Barnes-Hut tree. The hash below is a 
     simple interleaving of the bits of each component of a particle's
     position vector.

     The idea is to locate each particle within a three-dimensional grid which
     has BOXES_PER_DIM points in each dimension, and whose extents are
     the same as those of the simulation universe. Thus, for each particle, 
     we use its position coordinates to obtain three integer values telling which
     of those grid points it falls on. The width of each integer is BOXES_PER_DIM,
     and by interleaving the bits of the three integers, we can tell the position
     of the particle in the Barnes-Hut tree.
   */
  void hashParticleCoordinates(OrientedBox<double> universe){
    Key prepend;
    prepend = 1L;
    prepend <<= (TREE_KEY_BITS-1);

    /*
       Real xsz = universe.greater_corner.x-universe.lesser_corner.x;
       Real ysz = universe.greater_corner.y-universe.lesser_corner.y;
       Real zsz = universe.greater_corner.z-universe.lesser_corner.z;
     */
    Vector3D<Real> sz = universe.greater_corner-universe.lesser_corner;
    Vector3D<Real> prel;

    for(int i = 0; i < myNumParticles; i++){
      Particle p = myParticles[i];
      // Obtain the integer grid points on which the particle falls in each dimension

      prel = p.position-universe.lesser_corner;
      prel /= sz;
      prel *= (1.0*BOXES_PER_DIM);

      Key xint = (Key) prel.x;
      Key yint = (Key) prel.y;
      Key zint = (Key) prel.z;

      /*
         Key xint = ((Key) (((p->position.x-universe.lesser_corner.x)*(BOXES_PER_DIM*1.0))/xsz)); 
         Key yint = ((Key) (((p->position.y-universe.lesser_corner.y)*(BOXES_PER_DIM*1.0))/ysz)); 
         Key zint = ((Key) (((p->position.z-universe.lesser_corner.z)*(BOXES_PER_DIM*1.0))/zsz)); 
       */

      // Interleave bits
      Key mask = Key(0x1);
      Key k = Key(0x0);
      int shiftBy = 0;
      for(int j = 0; j < BITS_PER_DIM; j++){
        k |= ((zint & mask) <<  shiftBy);
        k |= ((yint & mask) << (shiftBy+1));
        k |= ((xint & mask) << (shiftBy+2));
        mask <<= 1;
        // minus 1 because mask itself has shifted
        // left by one position
        shiftBy += (NDIMS-1);
      }
      // Prepend the key with a '1' bit.
      k |= prepend; 
      p.key = k;
    }
  }


  /*
     Find out how many particles held by this PE
     fall under the root node (i.e. contribute the
     number of particles held by this PE)
   */
  void initHistogramParticles(){
    int rootDepth = 0;

    root = new Node<ForceData>(1,
                               rootDepth,
                               myParticles.get(0),
                               myNumParticles);
    /*
       The activeBins data structure keeps track of the
       active leaves (bins) iteration after iteration. 
       We begin with the root as the only active leaf. 
       When a worker is told by the master to split an
       active leaf, it partitions its particles among
       its children and contributes the number of particles
       held by each child to a reduction.
     */
    activeBins.addNewNode(root);

    // don't access myParticles through ckvec after this
    // anyway. these must be reset before this DM starts
    // to receive submitted particles from TPs placed on it
    myNumParticles = 0;
    myParticles.length() = 0;
  }

  void sendHistogram(){
    CkCallback cb = new CkCallback(CkReductionTarget(DataManager,receiveHistogram),0,this.thisgroup);
    Array<int> tmp = new Array<int>(activeBins.getNumCounts()]);
    for(int i = 0; i < activeBins.getNumCounts(); i++){
      tmp[i] = activeBins.getCount(i);
    }
    contribute(tmp,CkReduction.sum_int,cb);
    activeBins.reset();
    delete cb;
    delete tmp;
  }

  void doneDecomposition(){
    CkVec<Key> retractSites = new CkVec<Key>();
    int globalNumParticles = setGlobalParticleCounts(root);
    findRetractSites(root,retractSites);
    //CkPrintf("retractable nodes %d particles %d\n", retractSites.length(), globalNumParticles);
    /*
       Tell all PEs that there are no remaining active leaves,
       i.e. we have obtained a partitioning of particles on to
       tree pieces such that each tree piece gets no more than 
       a threshold (ppc) number of particles.
     */
    thisProxy@sendParticles(retractSites);
    // FIXME - do i have to delete this?
    delete retractSites;
  }

  // must set global num particles for each node before
  // we mark retractable nodes; otherwise, the retractable
  // nodes will not be obtained in DFS order
  int setGlobalParticleCounts(Node node){
    int globalParticleCount;
    if(node.getNumChildren() == 0){
      // must have already set global num particles for leaves
      globalParticleCount = node.getOwnerStart();
      CkAssert(globalParticleCount >= 0);
    }
    else{
      // if internal, process children first, so that
      // their global particle counts are available
      globalParticleCount = 0;
      globalParticleCount += setGlobalParticleCounts(node.getLeftChild());
      globalParticleCount += setGlobalParticleCounts(node.getRightChild());
      node.setOwnerStart(globalParticleCount);
    }

    return globalParticleCount;
  }

  void findRetractSites(Node node, CkVec<Key> sites){
    Real thresh = (DECOMP_TOLERANCE*(Real)(globalParams.ppc));
    // if this node has fewer particles than the threshold for splitting,
    // check whether its parent also sub-threshold. since we would like the
    // list of retractable nodes to be minimal, we should put a node in it
    // iff the node is retractable but its parent isn't
    // also, it is wasteful to mark leaves as retractable
    if(node.getNumChildren() > 0 && node.getOwnerStart() <= thresh){ 
      // FIXME - set parent correctly while refining
      Node parent = node.getParent();
      if(parent != null && 
          parent.getOwnerStart() > thresh) sites.push_back(node.getKey());
    }
    else if(node.getNumChildren() > 0){
      findRetractSites(node.getLeftChild(i),sites);
      findRetractSites(node.getRightChild(i),sites);
    }
  }

  void skipFlushParticles(){
    int pstart = 0;
    int pend = myNumParticles;
    for(int i = 1; i < numTreePieces; i++){
      Key check = treePieceKeys.get(i);
      int prev_start = pstart;
      pstart = binary_search_ge_v1(check,myParticles.getVec(),pstart,pend); 

      int prev_np = pstart-prev_start;

      sendParticleMsg(i-1, myParticles.get(prev_start), prev_np);
    }
    sendParticleMsg(numTreePieces-1,myParticles.get(pstart),pend-pstart);
  }

  /*
     As a result of the decomposition procedure, the
     particles held by each PE have been partitioned 
     among the tree pieces. It is now time to send these
     particles to their respective tree pieces.

     We do this by traversing the tree that
     was constructed during the decompostion by the 
     activeBins data structure.

     Returns number of tree pieces (leaves of the decomposition tree)
   */

  int flushParticles(FlushParticlesStruct fps){
    treePieceKeys.resize(0);
    return flushAndMark(root,0,fps);
  }

  int flushAndMark(Node node, int leafNum, FlushParticlesStruct fps){
    int
    if(fps.doRetract(node.getKey())){
      //CkPrintf("Retract node %llu\n", node->getKey());
      node.deleteBeneath();
      fps.advance();
    }

    node.setOwnerStart(leafNum);
    if(node.getNumChildren() > 0){
      leafNum = flushAndMark(node.getLeftChild(),leafNum,fps);
      // returned value is used by next child
      leafNum = flushAndMark(node.getRightChild(),leafNum,fps);
    }
    else{
      // sendParticlesToTreePiece whose root is 'node'
      sendParticlesToTreePiece(node,leafNum);
      leafNum++;
    }
    node.setOwnerEnd(leafNum);
    return leafNum;
  }

  /*
     Send the particles that you are holding for tree piece 'tp' to it.
   */
  void sendParticlesToTreePiece(Node nd, int tp) {
    CkAssert(nd.getNumChildren() == 0);
    int np = nd.getNumParticles();
    CkAssert(tp == treePieceKeys.length());
    treePieceKeys.push_back(Node<ForceData>.getParticleLevelKey(nd));
    //CkPrintf("(%d) DECOMP TP %d node %llu\n", CkMyPe(), tp, nd->getKey());

    sendParticleMsg(tp,nd.getParticles(),np);
  }

  void sendParticleMsg(int tp, Particle p, int np){
    ParticleArray arrayToSend = new ParticleArray(p,np);
    ParticleMsg msg = new ParticleMsg(arrayToSend);
    treePieceProxy[tp]@receiveParticles(msg);
    delete arrayToSend;
  }

  /*
     How many tree pieces are hosted on this PE? 
     This tells the PE how many tree pieces it should
     expect particle contributions from before starting
     to build the local tree.
   */
  void senseTreePieces(){
    localTreePieces.reset();
    CkLocMgr mgr = treePieceProxy@ckLocMgr();
    iterateOverLocMgr(mgr,localTreePieces);
  }

  /*
     Once the DM has collected particles from all the tree pieces hosted
     on this PE, it copies these particles into a contiguous buffer,
     sorts them and partitions them in-place in much the same way as 
     was done during domain decomposition. The result of this partitioning
     is a local tree whose leaves are buckets containing particles on
     the PE. The nodes of this tree are of six types:

     1. Bucket: leaf containing local particles
     2. EmptyBucket: an empty bucket
     3. Internal: a node whose every child is either Internal or Bucket or EmptyBucket

     4. RemoteBucket: leaf that is local to some other PE
     5. EmptyRemoteBucket: empty RemoteBucket
     6. Remote: a node whose entire subtree is non-local to this PE

     7. Boundary: a node that is none of the above, i.e. is "shared" between several
     PEs.
   */

  // Each node should have been marked with the 
  // min and max tree piece indices that it
  // hosts beneath it. Use this information to
  // obtain the roots of the 'localTreePieces' 
  // and build trees there. 
  // This procedure only goes
  // down to the depths of TreePiece roots; after that,
  // the singleBuildTree function is invoked.
  void buildTree(Node node, int pstart, int pend, int tpstart, int tpend){
    TB_DEBUG("(%d) pstart %d pend %d tpstart %d tpend %d node %lu\n", CkMyPe(), pstart, pend, tpstart, tpend, node.getKey());

    int np;
    nodeTable.insert(make_pair(node.getKey(),node));

    if(tpend <= tpstart){
      // No local tree piece under this node
      // It is Remote/RemoteBucket/RemoteEmptyBucket
      // Make request for remote node
      int numOwners = node.getOwnerEnd()-node.getOwnerStart();
      int requestOwner = node.getOwnerStart()+(rand()%numOwners);
      TB_DEBUG("(%d) requestMoments from tree piece %d for node %lu\n", CkMyPe(), requestOwner, node.getKey());

      CkEntryOptions opts = new CkEntryOptions();
      opts.setQueueing(CK_QUEUEING_IFIFO);
      opts.setPriority(REQUEST_MOMENTS_PRIORITY);
      treePieceProxy[requestOwner]@requestMoments(node.getKey(),CkMyPe(),opts);
      numMomentsRequested++;
      delete opts;

      // There are no particles on this PE under this node
      CkAssert(pstart == pend);
      node.setParticles(null,0);
      // Delete subtree beneath this node
      node.deleteBeneath();
      // Set type when requested moments are received
      // Don't tell parent that I'm done

      return;
    }
    else if(tpend-tpstart == 1 && (node.getOwnerEnd()-node.getOwnerStart()==1)){
      TreePieceDescriptor currentTP = localTreePieces.get(tpstart);
      np = currentTP.numParticles;
      TB_DEBUG("(%d) SINGLE LOCAL tree piece %d for node %lu\n", CkMyPe(), currentTP.index, node.getKey());
      // This is the first node that has 
      // a single tree piece beneath it.

      // Construct the entire tree underneath this node
      // and report the first particle index that doesn't
      // lie within it
      CkAssert(node.getOwnerStart() == currentTP.index);
      // There must be these many particles under the root of this TreePiece
      CkAssert(np == pend-pstart);

      if(np > 0) node.setParticles(myParticles.get(pstart),np);
      else node.setParticles(null,0);

      // Set the bucket indices for this TreePiece
      // If this is not the 0-th local TreePiece to have its bucket indices set,
      // it can use the last bucket (exclusive) of the one previous to it as its first bucket
      if(tpstart > 0) currentTP.bucketStartIdx = localTreePieces.get(tpstart-1).bucketEndIdx; 
      // To obtain the last bucket, set it to the first
      // bucket and increment each time a new bucket is encountered (in singleBuildTree)
      currentTP.bucketEndIdx = currentTP.bucketStartIdx;
      // Build the completely local tree underneath the root current TreePiece
      singleBuildTree(node,currentTP);
      // Since this node was completely local, its moments must have been computed
      node.setChildrenMomentsReady();
      // Tell parent node that one of its children has computed its moments
      notifyParentMomentsDone(node);

      // Set the root of the TreePiece
      currentTP.root = node;
      return;
    }

    np = pend-pstart;

    if(np > 0) node.setParticles(myParticles.get(pstart),pend-pstart);
    else node.setParticles(null,0);

    // If we are here, this node has at least one local tree piece underneath it
    CkAssert(tpend > tpstart);

    // Make sure that the range of tree pieces
    // is contained within this node; otherwise,
    // we shouldn't have made this call at all.
    CkAssert(localTreePieces.get(tpstart).index >= node.getOwnerStart());

    // Find the first local TreePiece that has
    // an index beyond the left child's last
    // contained TreePiece, i.e. equal to or after
    // the right child's first contained TreePiece 
    int firstTreePiece = node.getRightChild().getOwnerStart();
    Key particleTestKey = Node::getParticleLevelKey(node.getRightChild());
    int tp = binary_search_ge_v2(firstTreePiece,localTreePieces.get(0),tpstart,tpend);
    int pp = binary_search_ge_v1(particleTestKey,myParticles.get(0),pstart,pend);
    buildTree(node.getLeftChild(),pstart,pp,tpstart,pp);

    tpstart = tp;
    pstart = pp;

    buildTree(node.getRightChild(),pstart,pend,tpstart,tpend);

    if(node.allChildrenMomentsReady()){
      // All descendants were able to construct their subtrees
      // from data local to this PE: they must all be internal,
      // and so must this node. Note that this node cannot be
      // a Bucket or EmptyBucket, since it has more than one
      // local TreePiece underneath it
      node.setType(Internal);
      // Create node's moments from those of children
      node.getMomentsFromChildren();
      // Tell node's parent it is done building subtree
      notifyParentMomentsDone(node);
    }
    else{
      // Some descendants of node were unable to
      // construct their moments from PE-local data,
      // i.e. they are Remote: therefore, this node must be 
      // Boundary. It CANNOT be Remote since tpstart < tpend for it
      node.setType(Boundary);
    }
  }

  void notifyParentMomentsDone(Node node){
    // Respond to any requestors for the moments of 'node'
    TB_DEBUG("(%d) Notify parent flushing moments for node %lu\n", CkMyPe(), node.getKey());
    momentsResponseHelper(node);

    Node parent = node.getParent();
    if(parent != null) parent.childMomentsReady();
    else{
      // The root's moments have been computed,
      // i.e. all particles were internal to this PE.
      CkAssert(node.getKey() == Key(1));
      treeReady();
    }
  }

  // Recursive method to build the tree under a single TreePiece.
  // given an array of particles
  void singleBuildTree(Node node, TreePieceDescriptor tp){
    node.setOwners(tp.index,tp.index+1);
    //nodeTable[node->getKey()] = node;

    int np = node.getNumParticles();
    if(np <= ((Real)globalParams.ppb*BUCKET_TOLERANCE)){
      if(np == 0) node.setType(EmptyBucket);
      else node.setType(Bucket);

      node.getMomentsFromParticles();
      // Add to the list of buckets on this PE
      myBuckets.push_back(node);
      // Record the fact that current tree piece has another bucket
      tp.bucketEndIdx++;

      // since we might be reusing the tree, it could happen
      // that this node was internal (not a leaf) in the previous 
      // iteration. delete descendants, if there are any
      node.deleteBeneath();
    }
    else{
      node.setType(Internal);
      // Shouldn't not have allocated any children yet.
      // !doSkipDecomposition => node->getNumChildren() == 0
      CkAssert(doSkipDecomposition || (node.getNumChildren() == 0));
      // Partition the particles under this
      // node among the two children. To do this,
      // find the first particle that does not belong
      // underneath the left child, which is the same
      // as the first child that does belong under the
      // right child.
      if(node.getNumChildren() == 0) node.refine();
      else node.reuseRefine();

      Node child = node.getLeftChild();
      singleBuildTree(child,tp);
      nodeTable.insert(make_pair(child.getKey(),child));

      child = node.getRightChild();
      singleBuildTree(child,tp);
      nodeTable.insert(make_pair(child.getKey(),child));

      node.getMomentsFromChildren();
    }
  }

  Node lookupNode(Key k){
    return nodeTable.find(k);
  }

  void requestMoments(Key k, int replyTo){
    //pendingMoments[k].push_back(replyTo);
    CkVec<int> requests = pendingMoments.get(k);
    requests.push_back(replyTo);
    TB_DEBUG("(%d) received requestMoments from pe %d node %lu\n", CkMyPe(), replyTo, k);

    if(doneTreeBuild){
      Node node = lookupNode(k);
      CkAssert(node != null);

      boolean ready = node.allChildrenMomentsReady();
      TB_DEBUG("(%d) ready %d node %lu\n", CkMyPe(), ready, k);

      if(ready) momentsResponseHelper(node);
    }
  }

  void respondToMomentsRequest(Node node, CkVec<int> replyTo){
    MomentsExchangeStruct moments = new MomentsExchangeStruct(node);
    CkEntryOptions opts = new CkEntryOptions();
    opts.setQueueing(CK_QUEUEING_IFIFO);
    opts.setPriority(RECV_MOMENTS_PRIORITY);

    for(int i = 0; i < replyTo.length(); i++){
      TB_DEBUG("(%d) responding to %d with node %lu\n", CkMyPe(), ckVecRead(replyTo,i), node.getKey());
      thisProxy[ckVecRead(replyTo,i)]@receiveMoments(moments, opts);
    }
    replyTo.length() = 0;
    delete opts;
    delete moments;
  }

  void DataManager::receiveMoments(MomentsExchangeStruct moments){
    Node node = lookupNode(moments.key);
    CkAssert(node != null);
    CkAssert(!treeMomentsReady);

    numMomentsReceived++;
    // update moments of leaf and pass these on 
    // to parent recursively; if there are requests
    // for these nodes, respond to them
    updateLeafMoments(node,moments);
  }

  void updateLeafMoments(Node node, MomentsExchangeStruct data){
    copyMomentsToNode(node,data);
    TB_DEBUG("(%d) updateLeafMoments %lu\n", CkMyPe(), node.getKey());
    passMomentsUpward(node);
  }

  void momentsResponseHelper(Node node){
    CkVec<int> replyList = pendingMoments.find(node.getKey());
    if(replyList != null){
      respondToMomentsRequest(node,replyList);
      pendingMoments.erase(node.getKey());
    }
  }

  void passMomentsUpward(Node node){
    TB_DEBUG("(%d) passUp %lu\n", CkMyPe(), node.getKey());
    momentsResponseHelper(node);

    // FIXME - remember to check out setting of parent
    Node parent = node.getParent();
    if(parent == null){
      CkAssert(node.getKey() == 1);
      treeReady();
    }else{
      parent.childMomentsReady();
      TB_DEBUG("[%d] passup children ready %d for node %lu\n", CkMyPe(), parent.getNumChildrenMomentsReady(), parent.getKey());
      if(parent.allChildrenMomentsReady()){
        parent.getMomentsFromChildren();
        passMomentsUpward(parent);
      }
    }
  }

  void treeReady(){
    nodeLevelMerge();
  }

  void nodeLevelMerge(){
    if(globalParams.doPrintTree){
      string name("unmerged");
      doPrintTree(name);
    }
    // FIXME - do i have to delete this?
    treeMergerProxy@ckLocalBranch().submitPeRoot(new TreeRootContainer(root,CkMyPe()));
  }

  // FIXME - start here
  entry void doneNodeLevelMerge(TreeRootContainer container){
    root = container.root;
    if(globalParams.doPrintTree){
      string name("merged");
      if(CkMyPe()%numPesPerNode == 0) doPrintTree(name);
    }

    // register (i.e. fill into nodeTable)
    // nodes that are on paths to your useful 
    // tree pieces
    //CkPrintf("DM %d register\n", CkMyPe());
    registerTopLevelNodes(root,0,numLocalUsefulTreePieces);
    treeMomentsReady = true;
    flushBufferedRemoteDataRequests();

    CkCallback cb(CkReductionTarget(DataManager,doneTreeBuilds),0,thisProxy);
    contribute(cb);
  }

  void doneTreeBuilds(){
    double t = CmiWallTimer();
    CkPrintf("[0] Tree build took %f s\n", t-phaseTime);
    phaseTime = t;

    thisProxy@resumeDoneNodeLevelMerge();
  }

  void resumeDoneNodeLevelMerge(){
    startTraversal();
  }

  void registerTopLevelNodes(Node node, int tpstart, int tpend){
    nodeTable.insert(make_pair(node.getKey(),node));
    //nodeTable[node->getKey()] = node;
    if(tpend <= tpstart){
      //CkPrintf("%d REGISTER %llu external\n", CkMyPe(), node->getKey());
      return;
    }
    else if(tpend-tpstart == 1 && (node.getOwnerEnd()-node.getOwnerStart()==1)){
      TreePieceDescriptor tp = localTreePieces.get(tpstart);

      // if this node was a bucket and is being registered here,
      // it could be the merged replacement for a previous node.
      // therefore we must correct the pointer to this bucket in
      // myBuckets
      if(node.getType() == Bucket || node.getType() == EmptyBucket){
        // since this is the root of a tree piece (because registerTopLevelNodes
        // does not descend lower than the roots of tree pieces) and a bucket
        // it must be the only bucket assigned to the tree piece
        int bucketIdx = tp.bucketStartIdx;
        CkAssert(tp.bucketStartIdx+1 == tp.bucketEndIdx);
        //CkPrintf("%d REGISTER %llu bucket idx %d\n", CkMyPe(), node->getKey(), bucketIdx);
        ckVecWrite(myBuckets,bucketIdx,node);
      }

      tp.root = node;
    }
    else if(node.getNumChildren() > 0){
      Node leftChild = node.getLeftChild();
      Node rightChild = node.getRightChild();

      int firstGlobalTpUnderRight = rightChild.getOwnerStart();
      int idxFirstLocalTpUnderRight = binary_search_ge_v2(firstGlobalTpUnderRight,localTreePieces.get(0),tpstart,tpend); 

      registerTopLevelNodes(leftChild,tpstart,idxFirstLocalTpUnderRight);
      registerTopLevelNodes(rightChild,idxFirstLocalTpUnderRight,tpend);
    }
  }

  void flushBufferedRemoteDataRequests(){
    CkAssert(treeMomentsReady);
    for(int i = 0; i < bufferedNodeRequests.length(); i++){
      NodeRequest req = bufferedNodeRequests.get(i);
      requestNode(req.key, req.replyTo);
    }

    for(int i = 0; i < bufferedParticleRequests.length(); i++){
      NodeRequest req = bufferedParticleRequests.get(i);
      requestParticles(req.key, req.replyTo);
    }

    bufferedNodeRequests.length() = 0;
    bufferedParticleRequests.length() = 0;
  }

  void startTraversal(){
    LBTurnInstrumentOn();

    if(numLocalUsefulTreePieces > 0){
      int dummy=0;
      for(int i = 0; i < numLocalUsefulTreePieces; i++){
        TreePieceDescriptor tp = localTreePieces.get(i);
        tp.prepare(root,myBuckets);

        CkEntryOptions opts = new CkEntryOptions();
        opts.setPriority(START_TRAVERSAL_PRIORITY);
        opts.setQueueing(CK_QUEUEING_IFIFO);
        treePieceProxy[tp.index]@startTraversal(dummy,opts);
      }
    }
    else{
      finishIteration();
    }
  }

  ExternalParticle requestParticles(Node leaf, CutoffWorker worker, State state, Traversal traversal){
    Key key = leaf.getKey();
    ParticleCacheLine cacheLine = particleCache.get(key);
    if(cacheLine == null) cacheLine = particleCache.put(key);

    traceUserEvent(REMOTE_PARTICLE_REQUEST);

    if(cacheLine.data != null){
      return cacheLine.data;
    }
    else if(!cacheLine.sent){
      partReqs.incrRequests();

      cacheLine.sent = true;
      cacheLine.parentCached = leaf.isCached();
      if(leaf.isCached()){
        cacheLine.parent = leaf;
      }
      else{
        Node unsharedParent = leaf.copy(0);
        CkAssert(!unsharedParent.isCached());
        cacheLine.parent = unsharedParent;
      }

      CkAssert(leaf.getOwnerStart()+1 == leaf.getOwnerEnd());
      int owner = leaf.getOwnerStart();
      CkEntryOptions opts;
      opts.setQueueing(CK_QUEUEING_IFIFO);
      opts.setPriority(REQUEST_PARTICLES_PRIORITY);
      treePieceProxy[owner]@requestParticles(key, CkMyPe(), opts);

      RRDEBUG("(%d) REQUEST particles %lu leafCached %d from tp %d\n", CkMyPe(), key, leaf.isCached(), owner);
    }
    cacheLine.requestors.push_back(make_Requestor(worker,state,traversal,worker.getContext()));
    partReqs.incrDeliveries();
    return null;
  }

  void requestParticles(Key k, int replyTo) {
    if(!treeMomentsReady){
      bufferedParticleRequests.push_back(make_NodeRequest(k, replyTo));
      return;
    }

    RRDEBUG("(%d) REPLY particles key %lu to %d\n", CkMyPe(), k, replyTo);

    //map<Key,Node<ForceData>*>::iterator it = nodeTable.find(k);
    Node bucket = readTable(nodeTable,k);
    CkAssert(bucket != null);

    //CkAssert(it != nodeTable.end());
    //Node<ForceData> *bucket = it->second;
    CkAssert(bucket.getType() == Bucket);

    ExternalParticleArray externalParticleArray = makeExternalParticleArray(bucket);
    // FIXME - how to set priority?
    ParticleReplyMsg pmsg = new ParticleReplyMsg(k,externalParticleArray);
    /*
    *(int *)CkPriorityPtr(pmsg) = RECV_PARTICLES_PRIORITY;
    CkSetQueueing(pmsg,CK_QUEUEING_IFIFO);
    */

    thisProxy[replyTo]@recvParticles(pmsg);
  }

  Node requestNode(Node leaf, CutoffWorker worker, State state, Traversal traversal){
    Key key = leaf.getKey();
    NodeCacheLine nodeCacheLine = nodeCache.get(key);
    if(nodeCacheLine == null) nodeCacheLine = nodeCache.put(key);

    traceUserEvent(REMOTE_NODE_REQUEST);

    if(nodeCacheLine.data != null){
      RRDEBUG("(%d) REQUEST node %lu HIT!\n", CkMyPe(), key);
      return nodeCacheLine.data;
    }
    else if(!nodeCacheLine.sent){
      nodeReqs.incrRequests();

      nodeCacheLine.sent = true;
      nodeCacheLine.parentCached = leaf.isCached();
      if(leaf.isCached()){
        request.parent = leaf;
      }
      else{
        // in order to keep cached tree separate from 
        // local tree, since we can share trees across 
        // processors on an SMP node
        // by default, request.parentCached should be false
        request.parent = leaf.copy(0);
      }

      int numOwners = leaf.getOwnerEnd()-leaf.getOwnerStart();
      int requestOwner = leaf.getOwnerStart()+(rand()%numOwners);
      RRDEBUG("(%d) REQUEST node %lu leafCached %d from tp %d\n", CkMyPe(), key, leaf.isCached(), requestOwner);

      CkEntryOptions opts = new CkEntryOptions();
      opts.setPriority(REQUEST_NODE_PRIORITY);
      opts.setQueueing(CK_QUEUEING_IFIFO);
      treePieceProxy[requestOwner]@requestNode(key, CkMyPe(), opts);
      delete opts;
    }
    nodeCacheLine.requestors.push_back(make_Requestor(worker,state,traversal,worker.getContext()));
    nodeReqs.incrDeliveries();
    return null;
  }

  void requestNode(Key k, int replyTo){
    if(!treeMomentsReady){
      bufferedNodeRequests.add(make_NodeRequest(k,replyTo));
      return;
    }

    processNodeRequest(k,replyTo);
  }

  void processNodeRequest(Key key, int replyTo){
    CkAssert(treeMomentsReady);
    RRDEBUG("(%d) REPLY node %lu to %d\n", CkMyPe(), key, replyTo);

    //map<Key,Node<ForceData>*>::iterator it = nodeTable.find(key);
    Node node = readTable(nodeTable,key);
    CkAssert(node != null);
    //CkAssert(it != nodeTable.end());
    //Node<ForceData> *node = it->second;

    CkAssert(node.getNumChildren() > 0);

    // FIXME - how to set priority?
    NodeReplyMsg nmsg = node.serialize(globalParams.cacheLineSize);

    //*(int *)CkPriorityPtr(nmsg) = RECV_NODE_PRIORITY;
    //CkSetQueueing(nmsg,CK_QUEUEING_IFIFO);

    thisProxy[replyTo]@recvNode(nmsg);
  }

  void recvParticles(ParticleReplyMsg msg){
    ParticleCacheLine line = particleCache.get(msg.key);
    CkAssert(line != null);

    RRDEBUG("(%d) RECVD particle REPLY for key %lu\n", CkMyPe(), msg.key);

    CkAssert(line.requestors.length() > 0);
    CkAssert(line.sent);
    CkAssert(line.msg == NULL);

    line.msg = msg;

    // attach particles to bucket in tree 
    Node leaf = req.parent;
    CkAssert(leaf != null);
    CkAssert(leaf.getType() == RemoteBucket);
    leaf.setParticles(msg.particles.get(0),msg.particles.size());

    partReqs.decrRequests();
    partReqs.decrDeliveries(line.requestors.length());
    line.deliver();
  }

  void recvNode(NodeReplyMsg msg){
    NodeCacheLine line = nodeCache.find(msg.key);
    CkAssert(line != null);

    //Request &req = it->second;
    CkAssert(line.requestors.length() > 0);
    CkAssert(line.sent);
    CkAssert(line.msg == null);

    line.msg = msg;

    // attach recvd subtree to appropriate point in local tree

    Node node = line.parent;
    CkAssert(node != NULL);

    RRDEBUG("(%d) RECVD node REPLY for key %llu parent %llu cached %d\n", CkMyPe(), msg.key, node.getKey(), node.isCached());

    node.setLeftChild(msg.roots[0]);
    node.setRightChild(msg.roots[1]);

    nodeReqs.decrRequests();
    nodeReqs.decrDeliveries(line.requestors.length());

    line.deliverNode();
  }

  void traversalsDone(long pnInter, long ppInter, long openCrit){
    numTreePiecesDoneTraversals++;
    numInteractions[0] += pnInter;
    numInteractions[1] += ppInter;
    numInteractions[2] += openCrit;
    if(numTreePiecesDoneTraversals == numLocalUsefulTreePieces){
      finishIteration();
    }
  }

  void finishIteration(){
    CkCallback cb(CkReductionTarget(DataManager,doneForces),0,thisProxy);
    contribute(cb);
  }

  void doneForces(){
    double t = CmiWallTimer();
    CkPrintf("[0] Traversal took %f s\n", t-phaseTime);
    phaseTime = t;

    thisProxy@resumeFinishIteration();
  }

  void resumeFinishIteration(){
    if(CkMyPe()%numPesPerNode == 0){
      string name("merged");
      //doPrintTree(name);
    }


    //CkPrintf("DM %d finishIteration\n", CkMyPe());
    // can't advance particles here, because other PEs 
    // might not have finished their traversals yet, 
    // and therefore might need my particles

    CkAssert(nodeReqs.test());
    CkAssert(partReqs.test());

    //DtReductionStruct dtred;

    CkCallback cb(CkReductionTarget(DataManager,advance),thisProxy);
    contribute(numInteractions,CkReduction::sum_long,cb);

  }

  void advance(Array<long> interactions){

    //DtReductionStruct *dtred = (DtReductionStruct *)(msg->getData());

    kickDriftKick();
    //kickDriftKick(myBox.box,myBox.energy);

    Real pad = 0.00001;
    myBox.expand(pad);
    myBox.numParticles = myNumParticles;

    Real pn = interactions[0]/1e6;
    Real pp = interactions[1]/1e6;
    Real oc = interactions[2]/1e6;

    if(CkMyPe() == 0){
      CkPrintf("[STATS] node inter %.3f part inter %.3f open crit %.3f\n", pn, pp, oc);
    }

    iteration++;

    if(iteration == 13){
      traceBegin();
    }
    else if(iteration == 17){
      traceEnd();
    }

    if(iteration % globalParams.decompPeriod != 0){ 
      doSkipDecomposition = true;
    }
    else{
      doSkipDecomposition = false;
    }


    if(iteration == globalParams.iterations){
      CkCallback cb = CkCallback(CkReductionTarget(Main,niceExit),mainProxy);
      contribute(cb);
    }
    else if(iteration % globalParams.balancePeriod == 0){
      if(CkMyPe() == 0) CkPrintf("(%d) INITIATE LB\n", CkMyPe()); 
      for(int i = 0; i < localTreePieces.submittedParticles.length(); i++){
        TreePiece tp = localTreePieces.submittedParticles[i].owner;
        tp.startlb();
      }
      // must do decomposition after a load balancing step
      doSkipDecomposition = false;
    }
    else{
      init();
      for(int i = 0; i < localTreePieces.length(); i++){
        TreePieceDescriptor tpDescriptor = localTreePieces.get(i);
        TreePiece tp = tpDescriptor.owner;
        tp.cleanup();
      }

      CkCallback cb = CkCallback(CkReductionTarget(DataManager,recvUnivBoundingBox),thisProxy);
      contributeBoundingBox(cb);
      //contribute(sizeof(BoundingBox),&myBox,BoundingBoxGrowReductionType,cb);
    }
  }

  // FIXME - whose responsibility is it to delete input references?
  void recvUnivBoundingBox(Array<double> bounds){
    BoundingBox univBB = new BoundingBox();
    for(int i = 0; i < 3; i++){
      univBB.box.lesser_corner[i] = data[i];
      univBB.box.greater_corner[i] = -data[i+3];
    }
    univBB.pe = data[6];
    univBB.ke = data[7];

    decompose(univBB);
  }

  void freeCachedData(){
    particleCache.free();
    nodeCache.free();
  }

  void quiescence(){
    CkPrintf("QUIESCENCE dm %d pieces done %d (%d) nodereq %d partreq %d\n",
        CkMyPe(),
        numTreePiecesDoneTraversals,
        localTreePieces.count,
        nodeReqs.test(),
        partReqs.test()
        );

    CkCallback cb(CkReductionTarget(Main,quiescenceExit),mainProxy);
    contribute(cb);
  }

  void freeTree(){
    // delete the trees underneath your tree pieces
    // this can be done in parallel by all PEs on node
    for(int i = 0; i < numLocalUsefulTreePieces; i++){
      TreePieceDescriptor tp = localTreePieces.get(i);
      Node treePieceRoot = tp.root;
      CkAssert(treePieceRoot != null);
      treePieceRoot.deleteBeneath();
    }

    CkAssert(root != null);
    treeMergerProxy@ckLocalBranch().freeMergedTree();
  }

  void reuseTree(){
    CkAbort("Don't reuse decomposition tree: needs to be fixed for NODE_LEVEL_MERGE");
    /*
    Node<ForceData> tproot;
    for(int i = 0; i < numLocalUsefulTreePieces; i++){
      tproot = localTreePieces.submittedParticles[i].root;
      CkAssert(tproot != NULL);
      tproot->reuseTree();
    }

    if(root != NULL){
#ifdef NODE_LEVEL_MERGE
      treeMergerProxy.ckLocalBranch()->reuseMergedTree();
#else
      root->reuseTree();
#endif
    }
      */
  }

  void kickDriftKick(){
    Real particleEnergy;
    Real particleKinetic;
    Real particlePotential;

    myBox.reset();

    if(globalParams.doPrintAccel && (iteration == globalParams.iterations-1)){
      for(int i = 0; i < myParticles.length(); i++){
        Particle p = myParticles.get(i);
        CkPrintf("[%d] part %d p %f %f %f v %f %f %f\n", CkMyPe(), 
                  p.order, 
                  p.position.x, 
                  p.position.y, 
                  p.position.z,
                  p.velocity.x, 
                  p.velocity.y, 
                  p.velocity.z);
      }
    }

    if(iteration > 0){
      for(int i = 0; i < myParticles.length(); i++){
        Particle p = myParticles.get(i);
        p.velocity += globalParams.dthf*p.acceleration;
      }
    }


    for(int i = 0; i < myParticles.length(); i++){
      Particle p = myParticles.get(i);
      particlePotential = p.mass*p.potential;
      particleKinetic = 0.5*p.mass*p.velocity.lengthSquared();
      myBox.pe += particlePotential;
      myBox.ke += particleKinetic;

      p.position += globalParams.dtime*p.velocity 
        + 0.5*globalParams.dtime*globalParams.dtime*p.acceleration;

      p.velocity += globalParams.dthf*p.acceleration;

      p.interMass += p.mass;
      p.interMass -= univBox.mass;
      if(p.interMass < 0.0) p.interMass = -p.interMass;
      if(p.interMass > tolerance && bad < 100){
        CkPrintf("[%d] particle with interMass %g should be 0.0\n", CkMyPe(), p.interMass);
        bad++;
      }
      p.interMass = 0.0;

      myBox.grow(p.position);
      myBox.mass += p.mass;

      p.acceleration.x = 0.0;
      p.acceleration.y = 0.0;
      p.acceleration.z = 0.0;
      p.potential = 0.0;
    }

    CkAssert(bad == 0);
  }

  void init(){
    LBTurnInstrumentOff();
    CkAssert(pendingMoments.empty());
    // safe to reset here, since all tree pieces 
    // must have finished iteration
    freeCachedData();

    decompIterations = 0;

    doneTreeBuild = false;
    myBuckets.length() = 0;
    treeMomentsReady = false;
    numTreePiecesDoneTraversals = 0;
    numTreePiecesDoneRemoteRequests = 0;

    if(doSkipDecomposition){
      reuseTree();
    }
    else{
      freeTree();
      nodeTable.clear();
    }

    CkAssert(activeBins.getNumCounts() == 0);

    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;

    numMomentsRequested = numMomentsReceived = 0;
  }

  void resumeFromLB(){
    /* we delay the freeing of data structures to this point
     * because we want the tree pieces to be able to use the
     * tree, and in particular their roots for load balancing
     * */
    init();
    CkCallback cb(CkReductionTarget(DataManager,recvUnivBoundingBox),thisProxy);
    contributeBoundingBox(cb);
  }

  void contributeBoundingBox(CkCallback cb){
    Array<Real> data = new Array<Real>([8]);
    if(myNumParticles > 0){
      for(int i = 0; i < 3; i++){
        data[i] = myBox.box.lesser_corner[i];
        data[i+3] = -myBox.box.greater_corner[i];
      }
    }
    else{
      for(int i = 0; i < 3; i++){
        data[i] = LARGEST;
        data[i+3] = -SMALLEST;
      }
    }
    data[6] = myBox.pe;
    data[7] = myBox.ke;
    contribute(data,CkReduction.max_int,cb);
  }

  void pup(PUPer p){
    p|numRankBits;
    p|prevIterationStart;

    p|myParticles;
    p|myNumParticles;
    p|myBox;

    p|numTreePieces;
    p|iteration;

    p|compareEnergy;
  }

};



/* 
  Particles are decomposed onto TreePieces. Each
  TreePiece traverses the global tree (possibly creating
  communication requests) to compute the forces on its
  particles.
*/
array[1D] TreePiece {
  CkVec<ParticleMsg> decompMsgsRecvd = new CkVec<ParticleMsg>();
  int myNumParticles = 0;

  int iteration = 0;

  CkVec<Node> peBuckets;
  int myStartBucketIdx;
  int myEndBucketIdx;

  Node root;
  Node myRoot;

/*
  map<Key,set<Key> > bucketKeys;
  */

  LocalState localTraversalState = new LocalState();
  RemoteState remoteTraversalState = new RemoteState();

  LocalTraversalWorker localTraversalWorker = new LocalTraversalWorker();
  RemoteTraversalWorker remoteTraversalWorker = new RemoteTraversalWorker();

  Traversal trav = new Traversal();

  DataManager myDM = null;

  CkGroupID orbLBProxy;
  int numLB;
  boolean haveOrbLB;

  int totalNumTraversals;

  entry TreePiece(){
    usesAtSync = CmiTrue;
    findOrbLB();
    init();
  }

  entry void receiveParticles(ParticleMsg msg){
    decompMsgsRecvd.push_back(msg);
    myNumParticles += msg.numParticles;
  }

  entry void requestMoments(Key k, int replyTo){
    myDM.requestMoments(k,replyTo);
  }

  entry void startTraversal(int dummy){
    trav.setDataManager(myDM);

    if(myEndBucketIdx-myStartBucketIdx == 0){
      doneRemoteRequests();
      finishIteration();
      return;
    }

    remoteTraversalState.reset(this,myEndBucketIdx-myStartBucketIdx,ckVecRead(peBuckets,myStartBucketIdx+0));
    remoteTraversalWorker.reset(this,remoteTraversalState,ckVecRead(peBuckets,myStartBucketIdx+0));

    RescheduleMsg msg = new (NUM_PRIORITY_BITS) RescheduleMsg;
    SetMsgPriorityQueueing(msg,REMOTE_GRAVITY_PRIORITY);
    thisProxy[thisIndex]@doRemoteGravity(msg);

    localTraversalState.reset(this,myEndBucketIdx-myStartBucketIdx,ckVecRead(peBuckets,myStartBucketIdx+0));
    localTraversalWorker.reset(this,localTraversalState,ckVecRead(peBuckets,myStartBucketIdx+0));

    msg = new (NUM_PRIORITY_BITS) RescheduleMsg;
    SetMsgPriorityQueueing(msg,LOCAL_GRAVITY_PRIORITY);
    thisProxy[thisIndex]@doLocalGravity(msg);
  }

  entry void doLocalGravity(RescheduleMsg msg){
    int i;
    int limit;
    limit = myEndBucketIdx-myStartBucketIdx;
    for(i = 0; 
        i < globalParams.yieldPeriod && 
                   localTraversalState.currentBucketIdx < limit; 
        i++){
      trav.topDownTraversal(root,localTraversalWorker,localTraversalState);
      localTraversalState.currentBucketIdx++;
      localTraversalState.currentBucketPtr = ckVecRead(peBuckets,myStartBucketIdx+localTraversalState.currentBucketIdx);
      localTraversalWorker.setContext(localTraversalState.currentBucketPtr);
    }

    localTraversalState.decrPending(i);
    //CkPrintf("[%d] localGravity cur %d pending %d\n",thisIndex, localTraversalState.current, localTraversalState.pending);
    if(localTraversalState.current < limit){
      CkAssert(!localTraversalState.complete());
      thisProxy[thisIndex]@doLocalGravity(msg);
    }
    else{
      delete msg;
      if(localTraversalState.complete()){
        //CkPrintf("tree piece %d traversal done local\n", thisIndex);
        traversalDone();
      }
    }
  }

  entry void TreePiece::doRemoteGravity(RescheduleMsg msg){
    int i;
    int limit;
    limit = myEndBucketIdx-myStartBucketIdx;
    for(i = 0; 
        i < globalParams.yieldPeriod &&
              remoteTraversalState.currentBucketIdx < limit;
        i++){
      trav.topDownTraversal(root,remoteTraversalWorker,remoteTraversalState);
      remoteTraversalState.currentBucketIdx++;
      remoteTraversalState.currentBucketPtr = ckVecRead(peBuckets,myStartBucketIdx+remoteTraversalState.currentBucketIdx);
      remoteTraversalWorker.setContext(remoteTraversalState.currentBucketPtr);
    }

    remoteTraversalState.decrPending(i);
    // CkPrintf("[%d] remoteGravity cur %d pending %d\n", thisIndex, remoteTraversalState.current, remoteTraversalState.pending);
    if(remoteTraversalState.current < limit){
      CkAssert(!remoteTraversalState.complete());
      thisProxy[thisIndex]@doRemoteGravity(msg);
    }
    else{
      delete msg;
      if(remoteTraversalState.complete()){
        doneRemoteRequests();
        //CkPrintf("tree piece %d traversal done remote\n", thisIndex);
        traversalDone();
      }
    }
  }

  entry void requestParticles(NodeRequest request){
    // forward the request to the DM, since TPs don't
    // really own particles or nodes
    myDM.requestParticles(request);
  }

  entry void requestNode(NodeRequest request){
    myDM.requestNode(request);
  }

  entry void quiescence(){
    CkPrintf("QUIESCENCE tree piece %d proc %d submitted %d numBuckets %d trav_done %d outstanding local %d remote %d\n", 
        thisIndex,
        CkMyPe(),
        myNumParticles,
        myEndBucketIdx-myStartBucketIdx,
        totalNumTraversals,
        localTraversalState.pending,
        remoteTraversalState.pending);

    CkCallback cb(CkReductionTarget(Main,quiescenceExit),mainProxy);
    contribute(cb);
  }

  entry void doAtSync(){
    AtSync();
  }

  void ResumeFromSync() {
    //if(thisIndex == 0) CkPrintf("tree piece %d ResumeFromSync\n", thisIndex);
    contribute(CkCallback(CkReductionTarget(DataManager,resumeFromLB()),dataManagerProxy));
  }

  void init(){
    decompMsgsRecvd.length() = 0;
    totalNumTraversals = 2;
    myDM = dataManagerProxy@ckLocalBranch();
    root = null;
    peBuckets = null;
    myStartBucketIdx = -1;
    myEndBucketIdx = -1;
    myNumParticles = 0;
    myRoot = null;
    localTraversalState.finishedIteration();
    remoteTraversalState.finishedIteration();
  }

  void checkTraversals(){
  }

  void clearBucketsDebug(){
  }

  void prepare(Node _root, Node _myRoot, CkVec<Node> _peBuckets, int firstBucketIdx, int lastBucketIdx){

  {
    root = _root;
    myRoot = _myRoot;
    CkAssert(myRoot.getNumParticles() == myNumParticles);
    peBuckets = _peBuckets;
    myStartBucketIdx = myFirstBucketIdx;
    myEndBucketIdx = lastBucketIdx; 
  }

};

include "Node.h";


nodegroup TreeMerger {
  int numPesPerNode = CkNumPes()/CkNumNodes();
  int numSyncd = 0;

  CkVec<TreeRootContainer> myDataManagers = new CkVec<TreeRootContainer>();
  Node mergedRoot = null;

  public void init(){
    myDataManagers.resize(0);
  }

  public void submitPeRoot(TreeRootContainer container){
    CmiLock(__nodelock);
    CkAssert(container.root != null);
    myDataManagers.push_back(container);
    if(myDataManagers.size() == numPesPerNode){
      //CkPrintf("Merger %d merging\n", CkMyNode());
      CkVec<Node> toMerge = new CkVec<Node>();
      for(int i = 0; i < myDataManagers.length(); i++){
        TreeRootContainer container = ckVecRead(myDataManagers,i);
        toMerge.push_back(container.root);
        // FIXME do i have to do this?
        delete container;
      }

      CkAssert(toMerge.size() > 0);
      mergedRoot = merge(toMerge);
      delete toMerge;

      for(int i = 0; i < myDataManagers.length(); i++){
        // leave it to data managers to delete their previous roots
        dataManagerProxy[myDataManagers[i].pe]@doneNodeLevelMerge(new TreeRootContainer(mergedRoot));
      }
      //CkPrintf("Merger %d merge DONE\n", CkMyNode());
      init();
    }
    CmiUnlock(__nodelock);
  }

  private Node merge(CkVec<Node> toMerge){
    // there is only one valid node: return that; everything
    // beneath this node will be found only in the corresponding
    // tree
    Node toReturn = ckVecRead(toMerge,0);
    if(toMerge.size() == 1){
      return toReturn;
    }

    Node toReturn = pickBestNodeFromMergeSet(toMerge);

    CkVec<Node> nextLevelToMerge = new CkVec<Node>();

    nextLevelToMerge.length() = 0;
    // merge left children of merge set
    for(int i = 0; i < toMerge.size(); i++){
      Node mergeSetMember = ckVecRead(toMerge,i);
      if(mergeSetMember.getNumChildren() > 0){
        nextLevelToMerge.push_back(mergeSetMember.getLeftChild());
      }
    }
    toReturn.setLeftChild(merge(nextLevelToMerge));

    nextLevelToMerge.length() = 0;
    // merge left children of merge set
    for(int i = 0; i < toMerge.size(); i++){
      Node mergeSetMember = ckVecRead(toMerge,i);
      if(mergeSetMember.getNumChildren() > 0){
        nextLevelToMerge.push_back(mergeSetMember.getRightChild());
      }
    }
    toReturn.setRightChild(merge(nextLevelToMerge));
    // FIXME - parent is not set correctly

    delete nextLevelToMerge;

    // before returning the merged node from this 
    // level, go through the merge set and delete
    // all nodes other than 'toReturn'
    for(int i = 0; i < toMerge.size(); i++){
      Node mergeSetMember = ckVecRead(toMerge,i);
      if(mergeSetMember != toReturn) delete mergeSetMember;
    }
    return toReturn;
  }

  private Node pickBestNodeFromMergeSet(CkVec<Node> toMerge){
    Node toReturn = ckVecRead(toMerge,0);
    for(int i = 1; i < toMerge.size(); i++){
      Node mergeSetMember = ckVecRead(toMerge,i);
      if(hasBetterScore(mergeSetMember,toReturn)) toReturn = mergeSetMember;
    }
    return toReturn;
  }

  private boolean hasBetterScore(Node a, Node b){
    CkAssert(!(a.isLocal() && b.isLocal()));
    int aScore = score(a);
    int bScore = score(b);
    return aScore > bScore;
  }

  private static int SCORE_LOCAL = 2;
  private static int SCORE_BOUNDARY = 1;
  private static int SCORE_REMOTE = 0;

  private int score(Node node){
    if(node.isLocal()) return SCORE_LOCAL;
    if(node.isBoundary()) return SCORE_BOUNDARY;
    if(node.isRemote()) return SCORE_REMOTE;
  }

  public void freeMergedTree(){
    CmiLock(__nodelock);
    numSyncd++;
    if(numSyncd == numPesPerNode){
      //CkPrintf("Merger %d freeing\n", CkMyNode());
      numSyncd = 0;
      mergedRoot.deleteBeneath();
      delete mergedRoot;
      mergedRoot = null;
      //CkPrintf("Merger %d free DONE\n", CkMyNode());
    }
    CmiUnlock(__nodelock);
  }

  public void reuseMergedTree(){
    CmiLock(__nodelock);
    numSyncd++;
    if(numSyncd == numPesPerNode){
      //CkPrintf("Merger %d reusing\n", CkMyNode());
      numSyncd = 0;
      mergedRoot.reuseTree();
      //CkPrintf("Merger %d reuse DONE\n", CkMyNode());
    }
    CmiUnlock(__nodelock);
  }

};
