#include "defines.h"
#include "TreeMerger.h"

extern CProxy_DataManager dataManagerProxy;

TreeMerger::TreeMerger(){
  numPesPerNode = CkNumPes()/CkNumNodes();
  numSyncd = 0;
  mergedRoot = NULL;
}

void TreeMerger::init(){
  myDataManagers.resize(0);
}

void TreeMerger::submit(int pe, Node<ForceData>* root){
  CmiLock(__nodelock);
  CkAssert(root != NULL);
  myDataManagers.push_back(make_pair(pe,root));
  if(myDataManagers.length() == numPesPerNode){
    //CkPrintf("Merger %d merging\n", CkMyNode());
    CkVec<MergeStruct> toMerge;
    for(int i = 0; i < myDataManagers.length(); i++){
      root = myDataManagers[i].second;
      pe = myDataManagers[i].first;
      if(root->getChildren() != NULL) toMerge.push_back(MergeStruct(root,root->getChildren(),pe)); 
    }

    // might have no useful tree pieces on this SMP node
    if(toMerge.length() > 0){
      root = merge(toMerge);
      if(root->getNumChildren() > 0){
        bool allChildrenInternal = true;
        int numParticles = 0;
        Node<ForceData> *child = root->getChildren();
        for(int i = 0; i < root->getNumChildren(); i++){
          numParticles += child->getNumParticles();
          if(!child->isInternal()) allChildrenInternal = false;
          child++;
        }
        root->setNumParticles(numParticles);
        if(allChildrenInternal) root->setType(Internal);
        else root->setType(Boundary);
      }
    }

    mergedRoot = root;
    for(int i = 0; i < myDataManagers.length(); i++){
      // leave it to data managers to delete their previous roots
      dataManagerProxy[myDataManagers[i].first].doneNodeLevelMerge(PointerContainer(root));
    }
    //CkPrintf("Merger %d merge DONE\n", CkMyNode());
    init();
  }
  CmiUnlock(__nodelock);
}

Node<ForceData> *TreeMerger::merge(CkVec<MergeStruct> &toMerge){
  CkAssert(toMerge.length() > 0);
  if(toMerge.length() == 1){
    int key = toMerge[0].parent->getKey();
    int pe = toMerge[0].pe;
    //CkPrintf("return pe %d key %d\n", pe, key);
    return toMerge[0].parent;
  }

  // find a target merge buffer for this level
  // we pick the first one with the most number of particles underneath
  int pick = 0;
  Node<ForceData> *node = toMerge[pick].parent;
  int maxCount = node->getNumParticles();

  int count;
  for(int i = 1; i < toMerge.length(); i++){
    node = toMerge[i].child;
    count = node->getNumParticles();

    // second part of clause says that if even though an EmptyBucket 
    // will not have more particles than its peers (which should all be
    // RemoteEmptyBuckets) we should still pick it
    if(count > maxCount){
      maxCount = count;
      pick = i;
    }
  }

  // this the buffer we have picked as the target for the merge
  // we will return its parent
  Node<ForceData> *thisLevelBuffer = toMerge[pick].child;
  Node<ForceData> *thisLevelParent = toMerge[pick].parent;

  // merge the children of nodes in toMerge
  // this is the vector that is handed to the recursive calls
  CkVec<MergeStruct> newToMerge;
  // begin with leftmost child
  for(int j = 0; j < BRANCH_FACTOR; j++){
    // get merge candidates for next level
    for(int i = 0; i < toMerge.length(); i++){
      Node<ForceData> *node = (toMerge[i].child+j);
      if(node->getType() == Boundary || node->isInternal()){
        newToMerge.push_back(MergeStruct(node,node->getChildren(),toMerge[i].pe));
      }
    }
    if(newToMerge.length() > 0){
      // merge next level's candidates
      Node<ForceData> *mergedBufferParent = merge(newToMerge);
      // copy the values of the parent of the merged buffer
      // into the target of the merge
      Node<ForceData> *thisLevelNode = thisLevelBuffer+j;
      if(thisLevelNode != mergedBufferParent){
        copyNode(thisLevelNode,mergedBufferParent);
      }
      // update number of particles
      if(thisLevelNode->getNumChildren() > 0){
        bool allChildrenInternal = true;
        int numParticles = 0;
        Node<ForceData> *child = thisLevelNode->getChildren();
        for(int i = 0; i < thisLevelNode->getNumChildren(); i++){
          numParticles += child->getNumParticles();
          if(!child->isInternal()) allChildrenInternal = false;
          child++;
        }
        thisLevelNode->setNumParticles(numParticles);
        if(allChildrenInternal) thisLevelNode->setType(Internal);
        else thisLevelNode->setType(Boundary);
      }
    }
    // don't need to do anything when none of the trees are able
    // to offer a parent node with children: thisLevelNode doesn't
    // have to be updated with the subtree of one of its peers

    // reuse the memory for this list in the merging of the next child
    newToMerge.resize(0);
  }

  // delete the buffers that were not selected as merge targets
  for(int i = 0; i < toMerge.length(); i++){
    if(toMerge[i].child != thisLevelBuffer) delete[] toMerge[i].child;
  }

  return thisLevelParent;
}

void TreeMerger::copyNode(Node<ForceData> *target, Node<ForceData> *source){
  *target = *source;
  CkAssert(target->getKey() == source->getKey());
  CkAssert(target->getOwnerStart() == source->getOwnerStart());
  CkAssert(target->getOwnerEnd() == source->getOwnerEnd());
}

void TreeMerger::freeMergedTree(){
  CmiLock(__nodelock);
  numSyncd++;
  if(numSyncd == numPesPerNode){
    //CkPrintf("Merger %d freeing\n", CkMyNode());
    numSyncd = 0;
    mergedRoot->deleteBeneath();
    delete mergedRoot;
    mergedRoot = NULL;
    //CkPrintf("Merger %d free DONE\n", CkMyNode());
  }
  CmiUnlock(__nodelock);
}

void TreeMerger::reuseMergedTree(){
  CmiLock(__nodelock);
  numSyncd++;
  if(numSyncd == numPesPerNode){
    //CkPrintf("Merger %d reusing\n", CkMyNode());
    numSyncd = 0;
    mergedRoot->reuseTree();
    //CkPrintf("Merger %d reuse DONE\n", CkMyNode());
  }
  CmiUnlock(__nodelock);
}
