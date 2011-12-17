#ifndef TREE_MERGER_H
#define TREE_MERGER_H

#include "MeshStreamer.h"
#include "barnes.decl.h"
#include "Node.h"
#include "Descriptor.h"

#include <utility>
using namespace std;

class TreeMerger : public CBase_TreeMerger {
  int numPesPerNode;
  int numSyncd;

  CkVec<pair<int,Node<ForceData>*> > myDataManagers;
  Node<ForceData> *mergedRoot;


  void init();

  void copyNode(Node<ForceData> *target, Node<ForceData> *source);
  Node<ForceData> *merge(CkVec<pair<Node<ForceData>*,Node<ForceData>*> > &toMerge);

  public:
  
  TreeMerger();

  void submit(int pe, Node<ForceData>*root);
  void freeMergedTree();
  void reuseMergedTree();
};

#endif
