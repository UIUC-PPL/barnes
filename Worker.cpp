#include "Worker.h"
#include "DataManager.h"
#include "util.h"
#include "TreePiece.h"
#include "gravity.h"

extern CProxy_TreePiece treePieceProxy;

#if 0
int ParticleFlushWorker::work(Node<NodeDescriptor> *node){
  if(node->getNumChildren() > 0) return 1;

  dataManager->sendParticlesToTreePiece(node,leafCnt); 
  // Since this is a leaf and we are using Oct
  // decomposition, all the particles under this node
  // belong to the leafCnt-th TreePiece
  node->setOwners(leafCnt,leafCnt);
  leafCnt++;
  return 0;
}

int MomentsWorker::work(Node<ForceData> *node){
  int numChildren = node->getNumChildren();
  Particle *particles = node->getParticles();
  Node<ForceData> *parent = node->getParent();
  int numParticles = node->getNumParticles();

  if(numChildren == 0){
    // this node has no children, i.e. is a leaf
    // on this PE.
    setLeafType(node);

    if(node->getType() == Bucket || node->getType() == EmptyBucket){
      node->getMomentsFromParticles();
      if(parent != NULL) parent->childMomentsReady();
    }
  }
  else if(numChildren > 0){
    setTypeFromChildren(node);
    // don't have moments of remote nodes yet
    if(node->getType() == Internal){ 
      node->getMomentsFromChildren();
      if(parent != NULL) parent->childMomentsReady();
    }
  }
  nodeTable[node->getKey()] = node;

  // check whether this node is the root of a tree piece
  // on this PE
  int ostart = node->getOwnerStart();
  int oend = node->getOwnerEnd();
  // is the root owned by a single tree piece?
  if(node->getKey() == Key(1) && node->getType() == Internal && ostart == oend){
    tpRoots[ostart] = node;
  }
  // this node cannot be the root of a local TP
  else if(ostart != oend){
    // can its children be roots of some TPs?
    Node<ForceData> *child = node->getChildren();
    for(int i = 0; i < node->getNumChildren(); i++){
      int childOwnerStart = child->getOwnerStart();
      int childOwnerEnd = child->getOwnerEnd();
      NodeType childType = child->getType();
      if(childOwnerStart == childOwnerEnd
         && (childType == Internal || childType == Bucket || childType == EmptyBucket)){
        tpRoots[childOwnerStart] = child;
      }
      child++;
    }
  }
  return 1;
}

void MomentsWorker::setLeafType(Node<ForceData> *leaf){
  int ownerStart = leaf->getOwnerStart();
  int ownerEnd = leaf->getOwnerEnd();

  beg:
  // is the current PE the one that owns it?
  if(peTreePieces[curTP].index == ownerStart) {
    if(leaf->getNumParticles() > 0) leaf->setType(Bucket);
    else leaf->setType(EmptyBucket);

    buckets.push_back(leaf);
    return;
  }
  else if(ownerStart < peTreePieces[curTP].index){
    // type could have been set to EmptyBucket during construction
    if(leaf->getType() == Invalid){
      leaf->setType(Remote); 
      int numOwners = leaf->getOwnerEnd()-leaf->getOwnerStart()+1;
      int requestOwner = leaf->getOwnerStart()+(rand()%numOwners);
      TB_DEBUG("(%d) requestMoments %lu from tree piece %d\n", CkMyPe(), leaf->getKey(), requestOwner);

      CkEntryOptions opts;
      opts.setQueueing(CK_QUEUEING_IFIFO);
      opts.setPriority(REQUEST_MOMENTS_PRIORITY);
      treePieceProxy[requestOwner].requestMoments(leaf->getKey(),CkMyPe(),&opts);
    }
    return;
  }
  else do {
    curTP++;
  } while(ownerStart > peTreePieces[curTP].index);
  // :)
  goto beg;
}


void MomentsWorker::setTypeFromChildren(Node<ForceData> *node){
  Node<ForceData> *child = node->getChildren();
  int numChildren = node->getNumChildren();
  bool isInternal = true;
  for(int i = 0; i < numChildren; i++){
    NodeType childType = child->getType();
    if(childType != Internal 
       && childType != Bucket 
       && childType != EmptyBucket){
      isInternal = false; 
      break;
    }
    child++;
  }

  if(isInternal){
    // all children were either buckets or internal
    node->setType(Internal);
  }else{
    // at least one child was not a bucket and not internal
    // note that not all children can be non-bucket and non-internal
    // since we wouldn't have split a node on this PE if we
    // didn't have at least some particles under it
    // therefore, the node can't be remote
    node->setType(Boundary);
  }
}
#endif

int TraversalWorker::work(Node<ForceData> *node){
  NodeType type = node->getType();
  state->nodeEncountered(currentBucket->getKey(),node);
  bool keep = getKeep(type);

  if(!keep){
    state->nodeDiscarded(currentBucket->getKey(),node);
    return 0;
  }

  // basic opening criterion
  bool open = openCriterionBucket(node,currentBucket);
  state->incrOpenCriterion();
  if(open){
    state->nodeOpened(currentBucket->getKey(),node);
    return 1;
  }

  if(repeat(node->getType())) return 0;

  state->beforeForces(currentBucket,node->getKey());
  int computed = nodeBucketForce(node,currentBucket);
  state->nodeComputed(currentBucket,node->getKey());
  state->incrPartNodeInteractions(currentBucket->getKey(),computed);
  return 0;
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
bool RemoteTraversalWorker::repeat(NodeType type){
  if(type == Boundary) return true;
  return false;
}

void TraversalWorker::beforeParticleForces(Key k){
  state->beforeForces(currentBucket,k);
}

void TraversalWorker::work(ExternalParticle *particle){
  int computed = partBucketForce(particle,currentBucket);
  state->incrPartPartInteractions(currentBucket->getKey(),computed);
}

void TraversalWorker::bucketDone(Key k){
  state->bucketComputed(currentBucket,k);
}

bool LocalTraversalWorker::getKeep(NodeType type){
  return keep[type];
}

bool RemoteTraversalWorker::getKeep(NodeType type){
  return keep[type];
}

void RemoteTraversalWorker::done(){
  ownerTreePiece->doneRemoteRequests();
  //CkPrintf("tree piece %d traversal done remote resume\n", ownerTreePiece->getIndex());
  ownerTreePiece->traversalDone();
}

const bool LocalTraversalWorker::keep[] = {false,true,true,false,true,false,false,false};
const bool RemoteTraversalWorker::keep[] = {false,false,false,false,true,true,true,false};

int TreeSizeWorker::work(Node<ForceData> *node){
  int depth = node->getDepth();
  if(depth+1 <= cutoffDepth){
    numNodes += node->getNumChildren();
    return 1;
  }
  return 0;
}

