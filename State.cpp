#include "Node.h"
#include "State.h"
#include "TreePiece.h"

extern string NodeTypeString[];

void State::nodeEncountered(Key bucketKey, Node<ForceData> *node){
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) bucket %lu node %lu 0 encountered %s type %s radius %f\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucketKey, node->getKey(), getDescription().c_str(), NodeTypeString[node->getType()].c_str(), node->data.moments.rsq);
#endif
}

void State::nodeOpened(Key bucketKey, Node<ForceData> *node){
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) bucket %lu node %lu 1 opened %s type %s \n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucketKey, node->getKey(), getDescription().c_str(), NodeTypeString[node->getType()].c_str());
#endif
}

void State::nodeDiscarded(Key bucketKey, Node<ForceData> *node){
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) bucket %lu node %lu 2 discarding %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucketKey, node->getKey(), getDescription().c_str());
#endif

  insert(bucketKey,node->getKey());
}

void State::beforeForces(Node<ForceData> *bucket, Key k){
#ifdef VERBOSE_TRAVERSAL_INTERACTION
  //CkPrintf("beforeforces bucket %lu key %lu\n", bucket->getKey(), k);
  savedAccelerations.resize(bucket->getNumParticles());
  for(int i = 0; i < bucket->getNumParticles(); i++){
    savedAccelerations[i] = (bucket->getParticles()+i)->acceleration;
  }
#endif
}

void State::nodeComputed(Node<ForceData> *bucket, Key nodeKey){
#ifdef VERBOSE_TRAVERSAL_INTERACTION
  ostringstream oss;
  for(int i = 0; i < bucket->getNumParticles(); i++){
    Particle *p = bucket->getParticles()+i;
    Vector3D<Real> deltaAcc = p->acceleration-savedAccelerations[i];
    oss << deltaAcc.x << " " << deltaAcc.y << " " << deltaAcc.z << " ; " ;
  }
  CkPrintf("(%d,%d) bucket %lu node %lu 3 computing acc: %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucket->getKey(), nodeKey, oss.str().c_str());
#endif
  insert(bucket->getKey(),nodeKey);
}

void State::bucketComputed(Node<ForceData> *bucket, Key k){
#ifdef VERBOSE_TRAVERSAL_INTERACTION
  ostringstream oss;
  for(int i = 0; i < bucket->getNumParticles(); i++){
    Particle *p = bucket->getParticles()+i;
    Vector3D<Real> deltaAcc = p->acceleration-savedAccelerations[i];
    oss << deltaAcc.x << " " << deltaAcc.y << " " << deltaAcc.z << " ; "; 
  }
  CkPrintf("(%d,%d) bucket %lu node %lu 3 computing acc: %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucket->getKey(), k, oss.str().c_str());
#endif
  insert(bucket->getKey(),k);
}


