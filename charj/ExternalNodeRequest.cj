class ExternalNodeRequest {
  NodeReplyMsg msg = null;
  boolean sent = false;
  // FIXME - check entire code for places where references
  // are being stored in CkVec<>s. If memory is allocated
  // in the process, we must free it
  CkVec<Requestor> requestors = new CkVec<Requestor>();

  Node<ForceData> parent = null;
  boolean parentCached;

  /*
   * When the subtree under a requested node is obtained on the PE, find the 
   * workers that requested it, and perform a traversal on this subtree
   * for each one. Note that the worker may be in the midst of a traversal for
   * a different bucket from the one that generated this request. Therefore,
   * we must save and restore its current bucket through get/setContext.
   */
  void deliver(){
    for(int i = 0; i < requestors.length(); i++){
      Requestor req = requestors[i];
      CutoffWorker<ForceData> worker = req.worker;

      // Save a pointer to the current bucket of the worker.
      Node<ForceData> saveContext = worker.getContext();
      // Replace it with the one that generated the request
      // for this subtree
      worker.setContext(req.context);

      // Recall the kind of traversal the worker was performing
      // when it generated the request... 
      Traversal<ForceData> traversal = req.traversal; 
      // And a pointer to its current book-keeping state
      State state = req.state;
      // For each child of the parent whose children were requested...
      for(int j = 0; j < msg.getNumSubtreeRoots(); j++){
        // perform a top-down traversal of this subtree.
        traversal.topDownTraversal(msg.getSubtreeRoot(j),worker,state);
      }

      // Restore the current bucket of the worker
      worker.setContext(saveContext);
      // Check whether this traversal has completed.
      state.decrPending();
      if(state.complete()) worker.done();
    }
    requestors.clear();
  }
};
