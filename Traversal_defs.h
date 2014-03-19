#ifndef __TRAVERSAL_DEFS_H__
#define __TRAVERSAL_DEFS_H__

#include "DataManager.h"

template<typename T>
void Traversal<T>::startSphTraversal(Node<T> *root, CutoffWorker<T> *worker, State* st){
  
  CkPrintf("Starting Sph traversal [%d].\n", CkMyPe());
  NearNeighborState* state = (NearNeighborState*)st;
  stack<Node<T>*> ptrStk;
  stack<Node<T>*> leafStk;
  stack<Node<T>*> prntStk;
  int ret;

  ptrStk.push(root); 

  Node<T> *tmp[BRANCH_FACTOR];
  int numChildrenToExpand;
  //reach to the leaf nodes on the local tree and put them into the leaf stack
  while(!ptrStk.empty()){
    Node<T> *t_parent = ptrStk.top();
    Node<T> *t_children = t_parent->getChildren();
    int t_numChildren = t_parent->getNumChildren();
    ptrStk.pop();

    numChildrenToExpand = 0;
    for(int i = 0; i < t_numChildren; i++){	
      Node<T> *t_node = t_children+i;
      NodeType type = t_node->getType();
      if(type == Bucket)
        leafStk.push(t_node);
      else
        ptrStk.push(t_node);
    }
  }
  int i=0;
  //process each leaf node from the leaf stack, bottom up walk starts here
  while(!leafStk.empty()){
    CkPrintf("bottomUpTraversal for the leaf node %d. [%d]\n", (i++), CkMyPe());
    
    //get the top leaf from the stack and process it
    Node<T> *cur_leaf = leafStk.top();
    int ret = bottomUpTraversal(cur_leaf->parent, worker, st);
    if(ret) 
      CkPrintf("bottomUpTraversal for the leaf node %d is completed successfully. [%d]\n", (i++), CkMyPe());
    else CkPrintf("bottomUpTraversal for the leaf node %d is unfinished and will continute when remote requests come back. [%d]\n", (i++), CkMyPe());
    
    leafStk.pop();
  }
}

template<typename T>
int Traversal<T>::bottomUpTraversal(Node<T> *node, CutoffWorker<T> *worker, State* state){
  int ret = worker->work(node);
  //if open the node, process the node and go down first
  if(ret){
    processSphLeaf(node, worker, state);
    if(!topDownTraversalSph(node, worker, state)) return 0; //traversal suspended
  }
  //get the parent of the node
  Node<T>* parent = node->getParent();

  //move up in the tree until reaching the root or ballwithinbounds
  while(parent != NULL){
    ret = worker->work(parent); //parent does not need to call processSphLeaf
    if(ret){
      //expand the parent and do a top down walk now.
      int ret = topDownTraversalSph(parent, worker, state);
      if(!ret) return 0; //traversal suspended
      //move one step up
      parent = parent->getParent();
    }
  }
  return 1; //traversal successfully finished
}

template<typename T>
int Traversal<T>::topDownTraversalSph(Node<T> *root, CutoffWorker<T> *worker, State *state){
  stack<Node<T>*> ptrStk;
  ptrStk.push(root); 

  Node<T> *tmp[BRANCH_FACTOR];
  int numChildrenToExpand;

  //NOTE: this is a dept first search, bfs won't work.
  while(!ptrStk.empty()){
    Node<T> *t_parent = ptrStk.top();
    Node<T> *t_children = t_parent->getChildren();
    int t_numChildren = t_parent->getNumChildren();
    ptrStk.pop();

    numChildrenToExpand = 0;
    if(t_numChildren == 0){
      int ret = processSphLeaf(t_parent,worker,state);
      return 0;
    }

    for(int i = 0; i < t_numChildren; i++){
      Node<T> *t_node = t_children+i;
      int ret = worker->work(t_node);
      if(ret > 0){
        tmp[numChildrenToExpand] = t_node;
        numChildrenToExpand++;
      }
    }
    for(int i = numChildrenToExpand-1; i >= 0; i--){
      ptrStk.push(tmp[i]);
    }
  }
  return 1;
}


/*
template<typename T>
void Traversal<T>::bottomUpTraversal(Node<T> *root, CutoffWorker<T> *worker, State* st){
  CkPrintf("Starting bottumup traversal [%d].\n", CkMyPe());
  NearNeighborState* state = (NearNeighborState*)st;
  stack<Node<T>*> ptrStk;
  stack<Node<T>*> leafStk;
  stack<Node<T>*> prntStk;

  ptrStk.push(root); 

  Node<T> *tmp[BRANCH_FACTOR];
  int numChildrenToExpand;

  while(!ptrStk.empty()){
    Node<T> *t_parent = ptrStk.top();
    Node<T> *t_children = t_parent->getChildren();
    int t_numChildren = t_parent->getNumChildren();
    ptrStk.pop();

    numChildrenToExpand = 0;
    for(int i = 0; i < t_numChildren; i++){	
      Node<T> *t_node = t_children+i;
      NodeType type = t_node->getType();
      if(type == Bucket)
        leafStk.push(t_node);
      else
        ptrStk.push(t_node);
    }
  }
  int i=0;
  //process each leaf node from the leaf stack
  while(!leafStk.empty()){
    CkPrintf("Processing leaf node %d. [%d]\n", (i++), CkMyPe());
    Node<T> *cur_leaf = leafStk.top();
	//for each particle
	Particle* particles = cur_leaf->getParticles();
	for(int i=0; i<cur_leaf->getNumParticles(); i++){
		((NearNeighborState*)state)->startParticle();
		//process the local bucket
    		processSphLeaf(cur_leaf, state, &particles[i]);
		//move up to the parent until root or until ball within bounds
		Node<T>* parent = cur_leaf->getParent();
		while(parent != NULL){
		    //if(ballWithinBounds(parent->getBoundingBox(), cur_leaf->getCenter(), state->getCurrentRadius())){
		    //    break;
		    //}
		    //if a node while moving up is instersect-balls
		    if(intersect(parent->getBoundingBox(), cur_leaf->getCenter(), state->getCurrentRadius())){
		    //expand that node - traverse down to the leafs (except your own node)
		    
		    }
		    parent = parent->getParent();

		}
	}
    leafStk.pop();
  }
  worker->done();
}
*/

template<typename T>
int Traversal<T>::processSphLeaf(Node<T> *leaf, CutoffWorker<T> *worker, State *state){
  NodeType type = leaf->getType();

  CkAssert(type != Internal);
  CkAssert(type != Boundary);

  if(type == EmptyBucket || type == RemoteEmptyBucket) return 1;

  if(type == Bucket){
    Particle *particles = leaf->getParticles();
    int np = leaf->getNumParticles();
    CkAssert(particles);
    CkAssert(np > 0);
   
    for(int i = 0; i < np; i++){
      worker->work(particles+i);
    }
    worker->bucketDone(leaf->getKey());
  }
  else if(type == RemoteBucket){
    ExternalParticle *particles = (ExternalParticle *)leaf->getParticles();
    int np = leaf->getNumParticles();
    if(particles == NULL){
      state->incrPending();
      dm->requestParticles(leaf,worker,state,this);
      return 0; //remote request, walk needs to be suspended
    }

    for(int i = 0; i < np; i++){
      worker->work(particles+i);
    }
    worker->bucketDone(leaf->getKey());
  }
  else if(type == Remote){
    state->incrPending();
    dm->requestNode(leaf,worker,state,this);
    return 0; //remote request, walk needs to be suspended
  }
  return 1; //node is successfully processed
}


template<typename T>
void Traversal<T>::preorderTraversal(Node<T> *root, CutoffWorker<T> *worker){
  stack<Node<T>*> stk;
  stk.push(root);

  while(!stk.empty()){
    Node<T> *t_node = stk.top();
    stk.pop();

    int ret = worker->work(t_node);
    if(ret > 0){
      Node<T> *t_children = t_node->getChildren();
      for(int i = t_node->getNumChildren()-1; i >= 0; i--){
        stk.push(t_children+i);
      }
    }
  }
}

template<typename T>
void Traversal<T>::topDownTraversal_local(Node<T> *root, CutoffWorker<T> *worker){
  stack<Node<T>*> ptrStk;
  int ret = worker->work(root);
  if(ret == 0) return;

  ptrStk.push(root); 

  Node<T> *tmp[BRANCH_FACTOR];
  int numChildrenToExpand;

  while(!ptrStk.empty()){
    Node<T> *t_parent = ptrStk.top();
    Node<T> *t_children = t_parent->getChildren();
    int t_numChildren = t_parent->getNumChildren();
    ptrStk.pop();

    numChildrenToExpand = 0;
    for(int i = 0; i < t_numChildren; i++){
      Node<T> *t_node = t_children+i;
      ret = worker->work(t_node);
      if(ret > 0 && t_node->getNumChildren() > 0){
        tmp[numChildrenToExpand] = t_node;
        numChildrenToExpand++;
      } 
    }
    for(int i = numChildrenToExpand-1; i >= 0; i--){
      ptrStk.push(tmp[i]);
    }
  }
}

template<typename T>
void Traversal<T>::topDownTraversal(Node<T> *root, CutoffWorker<T> *worker, State *state){
  stack<Node<T>*> ptrStk;
  int ret = worker->work(root);
  if(ret == 0) return;

  ptrStk.push(root); 

  Node<T> *tmp[BRANCH_FACTOR];
  int numChildrenToExpand;

  while(!ptrStk.empty()){
    Node<T> *t_parent = ptrStk.top();
    Node<T> *t_children = t_parent->getChildren();
    int t_numChildren = t_parent->getNumChildren();
    ptrStk.pop();

    numChildrenToExpand = 0;
    if(t_numChildren == 0){
      processLeaf(t_parent,worker,state);
      continue;
    }

    for(int i = 0; i < t_numChildren; i++){
      Node<T> *t_node = t_children+i;
      ret = worker->work(t_node);
      if(ret > 0){
        tmp[numChildrenToExpand] = t_node;
        numChildrenToExpand++;
      }
    }
    for(int i = numChildrenToExpand-1; i >= 0; i--){
      ptrStk.push(tmp[i]);
    }
  }
}

template<typename T>
void Traversal<T>::processLeaf(Node<T> *leaf, CutoffWorker<T> *worker, State *state){
  NodeType type = leaf->getType();

  CkAssert(type != Internal);
  CkAssert(type != Boundary);

  if(type == EmptyBucket || type == RemoteEmptyBucket) return;

  if(type == Bucket){
    Particle *particles = leaf->getParticles();
    int np = leaf->getNumParticles();
    CkAssert(particles);
    CkAssert(np > 0);
   
    for(int i = 0; i < np; i++){
      worker->work(particles+i);
    }
    worker->bucketDone(leaf->getKey());
  }
  else if(type == RemoteBucket){
    ExternalParticle *particles = (ExternalParticle *)leaf->getParticles();
    int np = leaf->getNumParticles();
    if(particles == NULL){
      state->incrPending();
      dm->requestParticles(leaf,worker,state,this);
      return;
    }

    for(int i = 0; i < np; i++){
      worker->work(particles+i);
    }
    worker->bucketDone(leaf->getKey());
  }
  else if(type == Remote){
    state->incrPending();
    dm->requestNode(leaf,worker,state,this);
  }
}

template<typename T>
void Traversal<T>::postorderTraversal(Node<T> *root, CutoffWorker<T> *worker){
  // value returned from work() is ignored
  // since it doesn't make sense in a postorder
  // traversal
  stack<Node<T>*> cstk;
  stack<Node<T>*> nstk;

  cstk.push(root);
  nstk.push(root);

  while(!cstk.empty()){
    Node<T> *c_node = cstk.top();
    Node<T> *n_node = nstk.top();

    // we have processed the children of
    // n_node, but not n_node itself
    if(c_node != n_node){
      worker->work(n_node);
      nstk.pop();
      continue;
    }

    // here, we have a node whose children
    // haven't yet been processed

    cstk.pop();
    // first check whether the node has children
    if(c_node->getNumChildren() > 0){
      // we must explore the children of this node
      for(int i = c_node->getNumChildren()-1; i >= 0; i--){
        Node<T> *child = c_node->getChild(i);
        int num = child->getNumChildren();
        cstk.push(child);
        nstk.push(child);
      }
    } else {
      // this node does not have any children
      // so it is ok to visit it
      worker->work(n_node);
      nstk.pop();
    }
  }

  while(!nstk.empty()){
    worker->work(nstk.top());
    nstk.pop();
  }
}

#endif // __TRAVERSAL_DEFS_H__

