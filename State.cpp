#include "Node.h"
#include "State.h"
#include "TreePiece.h"

void State::nodeEncountered(Key bucketKey, Node<ForceData> *node){
#ifdef TREE_PIECE_LOG
  logFile << "bucket " << bucketKey 
    << " node " << node->getKey() 
    << " type " << NodeTypeString[node->getType()] 
    << " radius " << node->data.moments.radius
    << " cm " << node->data.moments.cm.x 
    << "," << node->data.moments.cm.y 
    << "," << node->data.moments.cm.z
    << endl;
#endif

#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) %s bucket %lu encountered node %lu type %s radius %f\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), description.c_str(), bucketKey, node->getKey(), NodeTypeString[node->getType()].c_str(), node->data.moments.radius);
#endif
}

void State::nodeOpened(Key bucketKey, Node<ForceData> *node){
#ifdef TREE_PIECE_LOG
#endif

#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) %s bucket %lu opened node %lu type %s \n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), description.c_str(), bucketKey, node->getKey(), NodeTypeString[node->getType()].c_str());
#endif
}

void State::nodeDiscarded(Key bucketKey, Node<ForceData> *node){
#ifdef TREE_PIECE_LOG
  logFile << "bucket " << bucketKey << " discarded node" << node->getKey() << endl;
#endif

#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) %s bucket %lu discarding node %lu\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), description.c_str(), bucketKey, node->getKey());
#endif

  insert(bucketKey,node->getKey());
}

void State::nodeComputed(Node<ForceData> *bucket, Key nodeKey){
#ifdef TREE_PIECE_LOG
  logFile << "bucket " << bucket->getKey() << " node " << nodeKey << " COMPUTE " << endl;
#endif
#ifdef VERBOSE_TRAVERSAL
  string nanStr = "";
  /*
  int numParticles = bucket->getNumParticles();
  Particle *particles = bucket->getParticles();
  for(int i = 0; i < numParticles; i++){
    if(isnan(particles[i].acceleration.length())){
      nanStr = "NAN DETECTED!\n";
      break;
    }
  }
  */
  CkPrintf("(%d,%d) %s bucket %lu (%s) computing node %lu\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), description.c_str(), bucket->getKey(), nanStr.c_str(), nodeKey);
#endif
  insert(bucket->getKey(),nodeKey);
}

void State::bucketComputed(Node<ForceData> *bucket, Key k){
#ifdef TREE_PIECE_LOG
  logFile << "bucket " << bucket->getKey() << " PARTICLES " << k << endl;
#endif
#ifdef VERBOSE_TRAVERSAL
  string nanStr = "";
  /*
  int numParticles = bucket->getNumParticles();
  Particle *particles = bucket->getParticles();
  for(int i = 0; i < numParticles; i++){
    if(isnan(particles[i].acceleration.length())){
      nanStr = "NAN DETECTED!";
      break;
    }
  }
  */
  CkPrintf("(%d,%d) %s bucket %lu (%s) computing bucket %lu\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), description.c_str(), bucket->getKey(), nanStr.c_str(), k);
#endif
  insert(bucket->getKey(),k);
}


