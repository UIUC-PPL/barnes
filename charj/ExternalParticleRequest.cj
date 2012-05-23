class ExternalParticleRequest {
  ParticleReplyMsg msg = null;
  boolean sent = false;
  // FIXME - check entire code for places where references
  // are being stored in CkVec<>s. If memory is allocated
  // in the process, we must free it
  CkVec<Requestor> requestors = new CkVec<Requestor>();

  Node<ForceData> parent = null;
  boolean parentCached;

  /*
   * When remote particles are received, we must supply them to all workers
   * on this PE that requested them. These workers will then use them for
   * computations that couldn't previously perform due to missing data.
   */
  void deliver(){
    // For each requestor of this set of particles
    for(int i = 0; i < requestors.length(); i++){
      Requestor req = requestors[i];
      // Which worker requested these particles?
      CutoffWorker<ForceData> worker = req.worker;

      /*
       * Workers maintain the current local bucket for which they
       * are performing the traversal (either remote or local). 
       * Currently, the worker could be in the middle of a traversal
       * for some other local bucket. We save the pointer to this 
       * local bucket here, and replace it with the one which had caused
       * the worker to request these particles in the first place.
       */
      Node<ForceData> saveContext = worker.getContext();
      worker.setContext(req.context);

      // so that the worker may do some book-keeping, if required 
      worker.beforeParticleForces(parent.getKey());
      for(int j = 0; j < msg.particles.length(); j++){ 
        // supply each particle in the bucket to the worker
        worker.work(msg.particles[j]);
      }
      // Tell the worker that it will receive no more particles from this
      // bucket
      worker.bucketDone(parent.getKey());

      // Restore the current bucket of the worker, in case
      // it is in the middle of some other traversal.
      worker.setContext(saveContext);

      // Update the corresponding state to reflect the fact that a previously
      // outstanding remote data request is now complete.
      State state = req.state;
      state.decrPending();
      // If there are no more outstanding data requests for the 
      // traversal, we are done with it.
      if(state.complete()) worker.done();
    }

    // No more requestors of this bucket to keep track of 
    requestors.clear();
  }
};
