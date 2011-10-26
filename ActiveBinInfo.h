#ifndef __ACTIVE_BIN_INFO_H__
#define __ACTIVE_BIN_INFO_H__

/*
 * CharmBH - ActiveBinInfo
 * This class is used to create a partition tree over a set of particles.
 * It is used (in slightly varying forms) in both the decomposition and
 * tree building phases. In the decomposition phase, it is used to partition
 * the particles present on each processor among tree pieces. In the
 * tree building phase, it creates a local tree out of the particles on 
 * a PE.
 *
 * An ActiveBinInfo object maintains the following data structures:
 * 1. oldvec: a list of tree nodes that are under consideration, i.e. might have to
 *    be partitioned if they do not meet a certain criterion. In decomposition, for 
 *    example, an node represents the prospective root of a tree piece, so that all
 *    particles beneath that node are to be assigned to the corresponding tree 
 *    piece. If this number is more than a threshold, we would like to partition
 *    the node. However, a PE knows only the number of *its* particles that lie 
 *    under each node in the partition tree. In order to determine whether a node
 *    needs to be partitioned, all PEs perform a reduction, so that one of them 
 *    (the master) knows which nodes should be partitioned, and which have fewer
 *    than the threshold number of particles. Therefore, the nodes in oldvec are those
 *    that are being considered for partitioning, not the ones that must be partitioned.
 *
 * 2. newvec: a list of newly obtained active nodes. These are the children of nodes that 
 *    have just been partitioned. This list of new active nodes is obtained by iterating through
 *    the oldvec list. For each node in the oldvec list, the master collects the global 
 *    number of particles beneath it. If this number is above a (user-specified, "ppc") threshold,
 *    the node is partitioned and its children are put into newvec.
 *
 * 3. counts: notice that the oldvec and newvec lists only maintain pointers to nodes; the actual
 *    count of number of particles beneath a node is held in the counts array. The counts array
 *    is a list of NodeDescriptor structs, which record not just the number of particles beneath 
 *    a node, but also the smallest and largest particle keys that have been found to exist within it.
 *
 * ActiveBinInfo is used in the following manner. The operations below are
 * performed in sequence on each PE, except for those marked with (!); these
 * are performed only by the master PE.
 *
 * 1. The root node is pushed into newvec through the addNewNode() method
 * 2. This also causes the number of particles beneath the root on this PE
 *    to be pushed into counts
 * 3. Use counts to contribute the total number of particles that this PE
 *    holds under each active node. The target of this reduction is the master
 *    PE. (see DataManager::sendHistogram)
 * 4. Reset, causing oldvec and newvec to be swapped
 * 5. (!) On the master, check which nodes need partitioning, and broadcast
 *    these decisions to all PEs (see DataManager::receiveHistogram).
 *    If none of the nodes need to be partitioned, inform all PEs that 
 *    decomposition has finished (see DataManager::sendParticles).
 * 6. Receive partitioning decisions from master, and apply them to the 
 *    nodes in the oldvec. (see DataManager::receiveSplitters)
 * 7. This causes the children of previously active nodes in oldvec
 *    to be pushed into the newvec list. Return to step 3.
 *    
 * 
 */

#include "charm++.h"
#include "Node.h"
#include <utility>
#include "Descriptor.h"

template<typename T>
struct ActiveBinInfo{ 
  CkVec<std::pair<Node<T>*, bool> > *oldvec;
  CkVec<std::pair<Node<T>*, bool> > *newvec;
  CkVec<int> counts;
  //CkVec<Node<T>*> unrefined;


  ActiveBinInfo(){
    oldvec = new CkVec<std::pair<Node<T> *, bool> >();
    newvec = new CkVec<std::pair<Node<T> *, bool> >();
  }

  ~ActiveBinInfo(){
    delete oldvec;
    delete newvec;
  }

  /*
   * Add a new node to be considered for partitioning. 
   * Used to inititate the decomposition/tree building 
   * process with the root of the tree.
   */
  void addNewNode(Node<T> *node){
    std::pair<Node<T>*,bool> pr;
    pr.first = node;
    pr.second = false;

    Key kfirst, klast;

    newvec->push_back(pr);
    int np = node->getNumParticles();
    Particle *particles = node->getParticles();
    /*
      Record the smallest and largest keys from
      particles on this PE that are under the root.
      Also put the number of such particles into
      counts.
    */
    if(np > 0){
      kfirst = particles[0].key;
      klast = particles[np-1].key;
    }else{
      kfirst = klast = Node<T>::getParticleLevelKey(node);
    }
    counts.push_back(np);
  }

  /*
   * This method is invoked to apply the partitioning decisions
   * on the nodes in oldvec. 
   */
  void processRefine(int *binsToRefine, int numBinsToRefine){
    // For each index specfied in the list of refinements
    for(int i = 0; i < numBinsToRefine; i++){
      // Get the current index of the node to refine
      int bin = binsToRefine[i];
      // Mark this node as refined.
      (*oldvec)[bin].second = true;
      Node<T> *node = (*oldvec)[bin].first;

      /*
       Actually refine the node; there are slightly different
       implementations of this method for decomposition and tree
       building, since tree building requires a little more work
       on each node (see OwnershipActiveBinInfo, below).
      */
      refine(node);

      Key kfirst, klast;

      std::pair<Node<T>*,bool> pr;
      pr.second = false;
      /*
       * For each child of the newly partitioned node, 
       * push the child into newvec so that it may be 
       * considered for partitioning next. 
       */
      for(int i = 0; i < node->getNumChildren(); i++){
        Node<T> *child = node->getChild(i);
        pr.first = child;
        
        // push child into newvec
        newvec->push_back(pr);
        int np = child->getNumParticles();
        Particle *particles = child->getParticles();
        // record smallest and largest keys on this
        // PE for each child
        if(np > 0){
          kfirst = particles[0].key;
          klast = particles[np-1].key;
        }else{
          kfirst = klast = Node<T>::getParticleLevelKey(child);
        }
        counts.push_back(np);
      }
    }

  }

  virtual void refine(Node<T> *node){
    node->refine();
  }

#if 0
  void processEmpty(int *emptyBins, int nEmptyBins){
    for(int i = 0; i < nEmptyBins; i++){
      int bin = emptyBins[i];
      Node<T> *node = (*oldvec)[bin].first;
      //node->core.type = NodeType::Empty;
    }
  }

  CkVec<Node<T>*> &getUnrefined(){
    for(int i = 0; i < oldvec->length(); i++){
      Node<T> *node = (*oldvec)[i].first;
      bool refined = (*oldvec)[i].second;

      if(!refined){
        unrefined.push_back(node);
      }
    }

    return unrefined;
  }
#endif

  /*
   * These are used when sending the histogram contributions
   * of this PE to the partitioning reduction.
   */
  int getNumCounts(){
    return counts.length();
  }

  int *getCounts(){
    return counts.getVec();
  }

  /*
   * Swap oldvec and newvec, and clear counts;
   * also make space for new children by clearing
   * newvec.
   */
  void reset(){
    CkVec<std::pair<Node<T>*, bool> > *tmp;
    tmp = oldvec;
    oldvec = newvec;
    newvec = tmp;
    newvec->length() = 0;

    //unrefined.length() = 0;
    counts.length() = 0;
  }

  CkVec<std::pair<Node<T>*, bool> > *getActive(){
    return oldvec;
  }
};


#if 0

bool CompareKeys(void *a, Key k);
/*
 * This is a customization of the ActiveBinInfo class for the
 * purpose of tree building. It differs from its parent in the
 * implementation of the refine() method. In addition to performing
 * basic refine() functions such as splitting the particles of a 
 * parent among its children, we also need to update the ownership
 * information of each node here, i.e. we need to track which tree
 * pieces own/share each node in the tree. This information is needed
 * when traversing the tree, so that if the particles on a PE need
 * to expand a remote node, they know which tree piece to ask for its
 * children.
 */
template<typename T>
struct OwnershipActiveBinInfo : public ActiveBinInfo<T> {

  /* 
   * This is a list of ranges of keys held by tree pieces.
   * As a result of the decomposition procedure, every tree
   * piece is assigned a certain range of particles. This list
   * is known to every PE (see keyRanges data field of DataManager)
   * The OwnershipActiveBinInfo class uses this information to
   * mark the owners of each node.
   *
   * The list is a series of [min_k,max_k] entries, where kmin and
   * kmax are the smallest and largest keys of particles held by
   * tree piece k.
   */
  Key *owners;

  OwnershipActiveBinInfo(Key *o) : 
    owners(o)
  {
  }

  void refine(Node<T> *node){
    /* 
     * See Node.h for definition of Node::refine()
     * This creates children for node, and initializes them with
     * the correct keys. It also partitions the particles held by
     * the parent among the children. 
     */
    node->refine();

    /*
     * In addition, we must mark the owner tree pieces of each 
     * node. The code below expects that the owners of the parent
     * have already been set. 
     */
    int ownerStart = node->getOwnerStart();
    int ownerEnd = node->getOwnerEnd();
    int start = (ownerStart<<1);
    int end = (ownerEnd<<1)+2;

    Node<T> *children = node->getChildren();
    Key childKey = children[0].getKey();

    /*
     * The first owner of the leftmost child is assumed to be the same as
     * that of the parent. This will be rectified if it turns out that the
     * first child is empty (see below, where owners are set to the garbage
     * value of -69,-69)
     */
    children[0].setOwnerStart(ownerStart);
    int childDepth = children[0].getDepth();
    childKey++;

    // For each child of the refined node,
    for(int i = 1; i < BRANCH_FACTOR; i++){
      // What is the smallest key of a particle that falls under this child?
      Key testKey = (childKey << (TREE_KEY_BITS-(childDepth*LOG_BRANCH_FACTOR+1)));
      /* 
       * In the keyRanges array, search for the first key that is GE the 
       * this smallest key. This gives us the range of owners for the child.
       */
      int owner_idx = binary_search_ge<Key>(testKey,owners,start,end,CompareKeys);
      // since each tree piece has two keys in the keyRanges
      // array (a lo and a hi), we divide by two to get the tp
      int tp_idx = (owner_idx>>1);

      // is the range of the tree piece contained completely
      // within the child node or does it straddle this child
      // and the previous one?
      if(EVEN(owner_idx)) children[i-1].setOwnerEnd(tp_idx-1);
      else                children[i-1].setOwnerEnd(tp_idx);

      if(children[i-1].getOwnerEnd() < children[i-1].getOwnerStart()){
        children[i-1].setOwners(-69,-69);
        children[i-1].setType(EmptyBucket);
      }

      children[i].setOwnerStart(tp_idx);
      start = owner_idx;
      childKey++;
    }

    // If it turns out that the last child is empty, set its owners to garbage
    children[BRANCH_FACTOR-1].setOwnerEnd(ownerEnd);
    if(children[BRANCH_FACTOR-1].getOwnerEnd() < children[BRANCH_FACTOR-1].getOwnerStart()){
      children[BRANCH_FACTOR-1].setOwners(-171,-171);
      children[BRANCH_FACTOR-1].setType(EmptyBucket);
    }
  }
};
#endif

#endif
