/*
 * CharmBH: Reductions.cpp
 * Customized reduction operations used in the program. 
 * The serial operations (grow(), operator+=, etc.) for each reduced 
 * data type (BoundingBox,NodeDescriptor, etc.) are defined in 
 * Descriptor.h. 
 */

#include "charm++.h"
#include "OrientedBox.h"
#include "defines.h"
#include "ActiveBinInfo.h"

CkReduction::reducerType BoundingBoxGrowReductionType;
#if 0
CkReduction::reducerType NodeDescriptorReductionType;
#endif
CkReduction::reducerType DtReductionType;

/*
 * Reduction type for constructing the bounding box of the entire
 * simulation universe. Each PE contributes the bounding box of its
 * particles.
 */
CkReductionMsg *BoundingBoxGrowReduction(int nmsgs, CkReductionMsg **msgs){
  Real *result = (Real *) (msgs[0]->getData());
  for(int i = 1; i < nmsgs; i++){
    CkAssert(msgs[i]->getSize() == 8*sizeof(Real));
    Real *compare = (Real *)msgs[i]->getData();
    for(int j = 0; j < 3; j++){
      if(compare[j] < result[j]) result[j] = compare[j];
      if(compare[j+3] > result[j+3]) result[j+3] = compare[j+3];
    }
    result[6] += compare[6];
    result[7] += compare[7];
  }

  return CkReductionMsg::buildNew(8*sizeof(Real),result);


  /*
  BoundingBox &bb = *((BoundingBox *) msgs[0]->getData());
  for(int i = 1; i < nmsgs; i++){
    CkAssert(msgs[i]->getSize() == sizeof(BoundingBox));
    BoundingBox &other = *((BoundingBox *) msgs[i]->getData());
    bb.grow(other);
  }
  return CkReductionMsg::buildNew(sizeof(BoundingBox),&bb);
  */
}

#if 0
/*
 * Reduction type for computing the total number of particles underneath 
 * active nodes during decomposition. Also computes the smallest and 
 * largest keys of particles under active nodes, across all PEs.
 * Each PE contributes the number of its particles that fall under an 
 * active node, as well as the smallest and largest keys among them.
 */
CkReductionMsg *NodeDescriptorReduction(int nmsgs, CkReductionMsg **msgs){
  int numElements = (msgs[0]->getSize()/sizeof(NodeDescriptor));

  for(int i = 0; i < numElements; i++){
    NodeDescriptor *zeroth = ((NodeDescriptor*)msgs[0]->getData())+i;
    for(int j = 1; j < nmsgs; j++){
      const NodeDescriptor &source = *(((NodeDescriptor *)msgs[j]->getData())+i);
      zeroth->grow(source);
    }
  }
  return CkReductionMsg::buildNew(msgs[0]->getSize(),msgs[0]->getData());
}
#endif

/*
 * Reduction type for basic statistics instrumented during the course of the
 * tree traversal on each PE. These are counts of particle-node interactions,
 * particle-particle interactions and number of calls to the opening criterion
 * function (gravity.h)
 */
CkReductionMsg *DtReduction(int nmsgs, CkReductionMsg **msgs){
  int numElements = (msgs[0]->getSize()/sizeof(DtReductionStruct));

  for(int i = 0; i < numElements; i++){
    DtReductionStruct *zeroth = ((DtReductionStruct*)(msgs[0]->getData()))+i;
    for(int j = 1; j < nmsgs; j++){
      const DtReductionStruct *source = ((const DtReductionStruct *)(msgs[j]->getData()))+i;
      *zeroth += *source;
    }
  }
  return CkReductionMsg::buildNew(msgs[0]->getSize(),msgs[0]->getData());
}

/*
 * Function called on each processor (initproc) to register the customized reduction
 * types with Charm++.
 */
void registerReducers(){
  BoundingBoxGrowReductionType = CkReduction::addReducer(BoundingBoxGrowReduction);
#if 0
  NodeDescriptorReductionType = CkReduction::addReducer(NodeDescriptorReduction);
#endif
  DtReductionType = CkReduction::addReducer(DtReduction);
}
