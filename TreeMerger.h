#ifndef TREE_MERGER_H
#define TREE_MERGER_H

#include "NDMeshStreamer.h"
#include "barnes.decl.h"
#include "Node.h"
#include "Descriptor.h"

#include <utility>
using namespace std;

struct MergeStruct {
  Node<ForceData> *parent;
  Node<ForceData> *child;
  int pe;

  MergeStruct(){
  }

  MergeStruct(Node<ForceData> *p, Node<ForceData> *c, int pe_) :
    parent(p), child(c), pe(pe_)
  {
  }
};

class TreeMerger : public CBase_TreeMerger {
  int numPesPerNode;
  int numSyncd;

  CkVec<pair<int,Node<ForceData>*> > myDataManagers;
  Node<ForceData> *mergedRoot;


  void init();

  void copyNode(Node<ForceData> *target, Node<ForceData> *source);
  Node<ForceData> *merge(CkVec<MergeStruct> &toMerge);

  public:
  
  TreeMerger();

  void submit(int pe, Node<ForceData>*root);
  void freeMergedTree();
  void reuseMergedTree();
};

#endif
