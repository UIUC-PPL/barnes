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

void State::nodeComputed(Node<ForceData> *bucket, Key nodeKey){
#ifdef VERBOSE_TRAVERSAL
  string nanStr = "";
  CkPrintf("(%d,%d) %s bucket %lu (%s) computing node %lu\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), getDescription().c_str(), bucket->getKey(), nanStr.c_str(), nodeKey);
#endif
  insert(bucket->getKey(),nodeKey);
}

void State::bucketComputed(Node<ForceData> *bucket, Key k){
#ifdef VERBOSE_TRAVERSAL
  string nanStr = "";
  CkPrintf("(%d,%d) %s bucket %lu (%s) computing bucket %lu\n", ownerTreePiece->getIndex(), ownerTreePiece->getIteration(), getDescription().c_str(), bucket->getKey(), nanStr.c_str(), k);
#endif
  insert(bucket->getKey(),k);
}


