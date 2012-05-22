class State {
  /*
   * Num. pending buckets; initially set to number of 
   * buckets owned by TreePiece, to indicate that traversals
   * need to be initiated for those buckets. Incremented by Traversal
   * whenever it encounters a remote node whose children
   * must be fetched, or a remote bucket whose particles must
   * be fetched. When nodes/particles are received, this
   * count is decreased. Also decremented when a TreePiece 
   * finishes a traversal with a bucket.
   */
  int pending;
  // Index of current bucket with which tree is being traversed
  // This field is used in TreePiece::doLocalGravity and 
  // TreePiece::doRemoteGravity 
  int currentBucketIdx;
  // Pointer to current bucket 
  Node<ForceData> currentBucketPtr;
  // The tree piece that initiated the traversal being tracked
  TreePiece *ownerTreePiece;

  // Number of interactions and opening tests performed during the traversal
  CmiUInt8 numInteractions[3];

  
  // Called when the traversal misses on remote data
  public void incrPending() { pending++; }

  // Either a remote request has been fulfilled or a bucket has finished traversal
  public void decrPending(int n=1){
    pending -= n;
  }

  // A traversal is complete when it has no pending remote data requests and
  // all buckets have finished traversal
  public bool complete() { return (pending == 0); }

  public State() : 
    pending(-1), 
    currentBucketIdx(-1), 
    currentBucketPtr(null),
    ownerTreePiece(null)
  {

    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;
  }
  
  // Reinitialize state
  public void reset(TreePiece owner, int p, Node<ForceData> bucketPtr){
  //void reset(TreePiece *owner, string &id, int p, Node<ForceData> **bucketPtr){
    ownerTreePiece = owner;
    //description = id;
    pending = p;
    currentBucketIdx = 0;
    currentBucketPtr = bucketPtr;
  }

  public void nodeEncountered(Key bucketKey, Node<ForceData> node){
    /*
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) bucket %lu node %lu 0 encountered %s type %s radius %f\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucketKey, node->getKey(), getDescription().c_str(), NodeTypeString[node->getType()].c_str(), node->data.moments.rsq);
#endif
  */
  }

  public void nodeOpened(Key bucketKey, Node<ForceData> node){
  /*
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) bucket %lu node %lu 1 opened %s type %s \n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucketKey, node->getKey(), getDescription().c_str(), NodeTypeString[node->getType()].c_str());
#endif
  */
  }

  public void nodeDiscarded(Key bucketKey, Node<ForceData> node){
    /*
#ifdef VERBOSE_TRAVERSAL
    CkPrintf("(%d,%d) bucket %lu node %lu 2 discarding %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucketKey, node->getKey(), getDescription().c_str());
#endif
    */

    //insert(bucketKey,node->getKey());
  }

  public void beforeForces(Node<ForceData> bucket, Key k){
    /*
#ifdef VERBOSE_TRAVERSAL_INTERACTION
    //CkPrintf("beforeforces bucket %lu key %lu\n", bucket->getKey(), k);
    savedAccelerations.resize(bucket->getNumParticles());
    for(int i = 0; i < bucket->getNumParticles(); i++){
      savedAccelerations[i] = (bucket->getParticles()+i)->acceleration;
    }
#endif
    */
  }

  void nodeComputed(Node<ForceData> bucket, Key nodeKey){
    /*
#ifdef VERBOSE_TRAVERSAL_INTERACTION
    ostringstream oss;
    for(int i = 0; i < bucket->getNumParticles(); i++){
      Particle *p = bucket->getParticles()+i;
      Vector3D<Real> deltaAcc = p->acceleration-savedAccelerations[i];
      oss << deltaAcc.x << " " << deltaAcc.y << " " << deltaAcc.z << " ; " ;
    }
    CkPrintf("(%d,%d) bucket %llu node %llu 3 computing acc: %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucket->getKey(), nodeKey, oss.str().c_str());
#endif
      */
    insert(bucket.getKey(),nodeKey);
  }

  void bucketComputed(Node<ForceData> bucket, Key k){
    /*
#ifdef VERBOSE_TRAVERSAL_INTERACTION
    ostringstream oss;
    for(int i = 0; i < bucket->getNumParticles(); i++){
      Particle *p = bucket->getParticles()+i;
      Vector3D<Real> deltaAcc = p->acceleration-savedAccelerations[i];
      oss << deltaAcc.x << " " << deltaAcc.y << " " << deltaAcc.z << " ; "; 
    }
    CkPrintf("(%d,%d) bucket %llu node %llu 3 computing acc: %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucket->getKey(), k, oss.str().c_str());
#endif
    */
    insert(bucket.getKey(),k);
  }


  // Clear counts of interactions/opening criterion applications
  public void finishedIteration(){
    numInteractions[0] = 0;
    numInteractions[1] = 0;
    numInteractions[2] = 0;
  }

    /*
  public void setBucketKeys(map<Key,set<Key> > bks){
#ifdef DEBUG_TRAVERSALS
    bucketKeys = bks;
#endif
  }
    */

  public void insert(Key bucketKey, Key k){
    /*
#ifdef DEBUG_TRAVERSALS
    //CkPrintf("%s interact bucket %llu node %llu\n", getDescription().c_str(), bucketKey, k);
    Key parentKey = (k>>1);
    set<Key>::iterator it = (*bucketKeys)[bucketKey].find(parentKey);
    while(it != (*bucketKeys)[bucketKey].end()){
      (*bucketKeys)[bucketKey].erase(it);
      parentKey >>= 1;
      it = (*bucketKeys)[bucketKey].find(parentKey);
    }
    (*bucketKeys)[bucketKey].insert(parentKey);
#endif
    */
  }

  // Methods to increment different statistics
  public void incrPartNodeInteractions(CmiUInt8 n){
    numInteractions[0] += n;
  }

  public void incrPartPartInteractions(CmiUInt8 n){
    numInteractions[1] += n;
  }

  public void incrOpenCriterion(){
    numInteractions[2]++;
  }

  public string getDescription() { return ""; }

  /*
#ifdef VERBOSE_TRAVERSAL_INTERACTION
  CkVec<Vector3D<Real> > savedAccelerations;
#endif
  */
};

class LocalState : State {
  string getDescription() { return "LOCAL"; }
};

class RemoteState : State {
  string getDescription() { return "REMOTE"; }
};


