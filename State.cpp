#include "Node.h"
#include "State.h"
#include "TreePiece.h"

extern string NodeTypeString[];

void State::nodeEncountered(Key bucketKey, Node<ForceData> *node){
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) %s bucket %lu encountered node %lu type %s radius %f\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), getDescription().c_str(), bucketKey, node->getKey(), NodeTypeString[node->getType()].c_str(), node->data.moments.rsq);
#endif
}

void State::nodeOpened(Key bucketKey, Node<ForceData> *node){
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) %s bucket %lu opened node %lu type %s \n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), getDescription().c_str(), bucketKey, node->getKey(), NodeTypeString[node->getType()].c_str());
#endif
}

void State::nodeDiscarded(Key bucketKey, Node<ForceData> *node){
#ifdef VERBOSE_TRAVERSAL
  CkPrintf("(%d,%d) %s bucket %lu discarding node %lu\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), getDescription().c_str(), bucketKey, node->getKey());
#endif

  insert(bucketKey,node->getKey());
}

void State::beforeForces(Node<ForceData> *bucket, Key k){
#if 0
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
  CkPrintf("(%d,%d) bucket %lu computing node %lu acc: %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucket->getKey(), nodeKey, oss.str().c_str());
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
  CkPrintf("(%d,%d) bucket %lu computing bucket %lu acc: %s\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), bucket->getKey(), k, oss.str().c_str());
#endif
  insert(bucket->getKey(),k);
}


