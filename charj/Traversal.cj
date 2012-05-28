class Traversal {

  // to request particles/nodes from
  DataManager dm = null;

  Node processLeaf(Node node, CutoffWorker worker, State state);

  public:

  void setDataManager(DataManager mgr){
    dm = mgr;
  }

  void topDownTraversal_local(Node root, CutoffWorker worker){
    if(worker.work(root)){
      for(int i = 0; i < root.getNumChildren(); i++){
        topDownTraversal_local(root.getChild(i),worker);
      }
    }
  }

  void topDownTraversal(Node root, CutoffWorker worker, State state){
    if(worker.work(root)){
      if(root.getNumChildren() == 0){
        // we might have a leaf that is a remote node, such that
        // the children of the remote node have been fetched for this
        // PE previously
        CkVec<Node> children = processLeaf(root,worker,state);
        if(children != null){
          for(int i = 0; i < children.size(); i++){
            topDownTraversal(children[i],worker,state);
          }
        }
      }
      else{
        for(int i = 0; i < root.getNumChildren(); i++){
          topDownTraversal(root.getChild(i),worker);
        }
      }
    }
  }

  Node<T> processLeaf(Node<T> leaf, CutoffWorker<T> worker, State state){
    NodeType type = leaf.getType();

    CkAssert(type != Internal);
    CkAssert(type != Boundary);

    if(type == EmptyBucket || type == RemoteEmptyBucket) return null;

    if(type == Bucket){
      CkVec<Particle> particles = leaf.getParticles();
      CkAssert(particles != null);
      CkAssert(particles.size() > 0);

      worker.beforeParticleForces(leaf.getKey());
      for(int i = 0; i < particles.size(); i++){
        worker.work(particles[i]);
      }
      worker.bucketDone(leaf.getKey());
      return null;
    }
    else if(type == RemoteBucket){
      CkAssert(leaf.getParticles() == null);
      CkVec<ExternalParticle> cached = dm.requestParticles(leaf,worker,state,this);
      if(cached != null){ 
        worker.beforeParticleForces(leaf.getKey());
        for(int i = 0; i < cached.size(); i++){
          worker.work(cached+i);
        }
        worker.bucketDone(leaf.getKey());
      }
      else{
        state.incrPending();
      }
      return null;
    }
    else if(type == Remote){
      CkVec<Node<T> > cached = dm.requestNode(leaf,worker,state,this);
      if(cached == null) state.incrPending();
      return cached;
    }

    CkAbort("bad leaf type\n");
    return null;
  }
};


