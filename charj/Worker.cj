// worker abstract class
abstract class CutoffWorker {
  public abstract int work(Node);
  public void work(ExternalParticle) {}
  public void bucketDone(Key k) {}
  public Node getCurrentBucket() {return null;}
  public void setCurrentBucket(Node bucket) {}
  public void done() {}
  public void beforeParticleForces(Key k) {}
};

class TraversalWorker : public CutoffWorker {
  TreePiece ownerTreePiece;
  State state;
  
  Node currentBucket;

  public TraversalWorker(){
    ownerTreePiece = null;
    currentBucket = null;
  }

  public void reset(TreePiece owner, State s, Node bucket){
    ownerTreePiece = owner;
    state = s;
    currentBucket = bucket;
  }

  public Node getCurrentBucket(){
    return currentBucket;
  }

  public void setCurrentBucket(Node current){
    currentBucket = current;
  }

  public int work(Node node){
    NodeType type = node.getType();
    state.nodeEncountered(currentBucket.getKey(),node);
    boolean keep = getKeep(type);

    if(node.getType() == EmptyBucket ||
        node.getType() == RemoteEmptyBucket){
      state.nodeComputed(currentBucket,node.getKey());
      return 0;
    }
    if(!keep){
      state.nodeDiscarded(currentBucket.getKey(),node);
      return 0;
    }

    // basic opening criterion
    boolean open;
    open = openCriterionBucket(node,currentBucket);
    state.incrOpenCriterion();
    if(open){
      state.nodeOpened(currentBucket.getKey(),node);
      return 1;
    }

    if(repeat(node.getType())) return 0;

    int computed;
    state.beforeForces(currentBucket,node.getKey());
    computed = nodeBucketForce(node,currentBucket);
    state.nodeComputed(currentBucket,node.getKey());
    state.incrPartNodeInteractions(computed);
    return 0;
  }

  public void work(ExternalParticle particle){
    int computed;
    computed = partBucketForce(particle,currentBucket);
    state.incrPartPartInteractions(computed);
  }

  public void bucketDone(Key k){
    state.bucketComputed(currentBucket,k);
  }
  
  public abstract boolean getKeep(NodeType type);
  public boolean repeat(NodeType type) {return false;}


  public void beforeParticleForces(Key k){
    state.beforeForces(currentBucket,k);
  }
};

class LocalTraversalWorker : TraversalWorker {
  private final static boolean keep[] = {false,true,true,false,true,false,false,false};

  public LocalTraversalWorker(){
    super();
  }

  public boolean getKeep(NodeType type){
    return keep[type];
  }
};

class RemoteTraversalWorker : public TraversalWorker {
  static final boolean keep[] = {false,false,false,false,true,true,true,false};
  
  public RemoteTraversalWorker(){
    super();
  }

  public void done(){
    ownerTreePiece.doneRemoteRequests();
    ownerTreePiece.traversalDone();
  }

  public boolean getKeep(NodeType type){
    return keep[type];
  }

/* Both remote and local workers process Boundary nodes 
 * However, if a Boundary node is far enough from a bucket
 * (if so, the repeat() method is called)
 * only a local worker should compute with it, 
 * not a remote worker. this avoids duplicate computation.
 * Also, the remote worker shouldn't skip any other
 * type of node that it processes, since the local
 * worker will not access it.
 * *
 * */
  public boolean repeat(NodeType type){
    if(type == Boundary) return true;
    return false;
  }
};

class TreeSizeWorker : CutoffWorker {
  int numNodes;
  int cutoffDepth;

  public TreeSizeWorker(){}

  public reset(int d){ 
    cutoffDepth = d;
    numNodes = 0;
  }

  public int work(Node node){
    int depth = node.getDepth();
    if(depth+1 <= cutoffDepth){
      numNodes += node.getNumChildren();
      return 1;
    }
    return 0;
  }

  public int getNumNodes(){
    return numNodes;
  }
};


