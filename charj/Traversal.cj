class Traversal {

  // to request particles/nodes from
  DataManager dm = null;

  Node processLeaf(Node node, CutoffWorker worker, State state);

  void setDataManager(DataManager mgr){
    dm = mgr;
  }

  void topDownTraversal_local(Node node, CutoffWorker worker){
    if(worker.work(node) && node.getNumChildren() > 0){
      topDownTraversal_local(node.getLeftChild(),worker);
      topDownTraversal_local(node.getRightChild(),worker);
    }
  }

  void topDownTraversal(Node node, CutoffWorker worker, State state){
    if(!worker.work(node)) return;

    if(node.getNumChildren() == 0){
      // we might have a leaf that is a remote node, such that
      // the children of the remote node have been fetched for this
      // PE previously
      node = processLeaf(node,worker,state);
      if(node == null) return;
      CkAssert(node.getNumChildren() > 0);
    }
    topDownTraversal(node.getLeftChild(),worker,state);
    topDownTraversal(node.getRightChild(),worker,state);
  }

  Node processLeaf(Node leaf, CutoffWorker worker, State state){
    NodeType type = leaf.getType();

    CkAssert(type != Internal);
    CkAssert(type != Boundary);

    if(type == EmptyBucket || type == RemoteEmptyBucket) return null;

    if(type == Bucket){
      CkAssert(leaf.getParticles() != null);
      CkAssert(leaf.getNumParticles() > 0);

      worker.beforeParticleForces(leaf.getKey());
      iterateOverParticles(leaf.getParticles(),leaf.getNumParticles(),worker);
      worker.bucketDone(leaf.getKey());
      return null;
    }
    else if(type == RemoteBucket){
      CkAssert(leaf.getParticles() == null);
      Node cachedLeaf = dm.requestParticles(leaf,worker,state,this);
      if(cachedLeaf != null){ 
        worker.beforeParticleForces(leaf.getKey());
        iterateOverExternalParticles(cachedLeaf.getParticles(), cachedLeaf.getNumParticles(), worker);
        worker.bucketDone(leaf.getKey());
      }
      else{
        state.incrPending();
      }
      return null;
    }
    else if(type == Remote){
      Node cachedParent = dm.requestNode(leaf,worker,state,this);
      if(cachedParent == null) state.incrPending();
      return cachedParent;
    }

    CkAbort("bad leaf type\n");
    return null;
  }
};


