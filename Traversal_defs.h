#ifndef __TRAVERSAL_DEFS_H__
#define __TRAVERSAL_DEFS_H__

#include "DataManager.h"
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
  //CkPrintf("(%d) VISIT node %ld\n", CkMyPe(), getKey());
  int ret = worker->work(root);
  if(ret == 0) return;

  ptrStk.push(root); 
  //CkPrintf("(%d) PUSH node %ld\n", CkMyPe(), getKey());

  Node<T> *tmp[BRANCH_FACTOR];
  int numChildrenToExpand;
  //tmp.reserve(2);

  while(!ptrStk.empty()){
    //pair<Node<T>*,int> pr = ptrStk.top();
    Node<T> *t_parent = ptrStk.top();
    Node<T> *t_children = t_parent->getChildren();
    int t_numChildren = t_parent->getNumChildren();
    //CkPrintf("(%d) POP node %ld\n", CkMyPe(), t_parent->getKey());
    ptrStk.pop();

    //tmp.length() = 0;
    numChildrenToExpand = 0;
    for(int i = 0; i < t_numChildren; i++){
      Node<T> *t_node = t_children+i;
      //CkPrintf("(%d) VISIT node %ld\n", CkMyPe(), t_node->getKey());
      ret = worker->work(t_node);
      if(ret > 0 && t_node->getNumChildren() > 0){
        //tmp.push_back(t_node);
        tmp[numChildrenToExpand] = t_node;
        numChildrenToExpand++;
      } 
    }
    for(int i = numChildrenToExpand-1; i >= 0; i--){
      //CkPrintf("(%d) PUSH child %d key %ld\n", CkMyPe(), i, tmp[i]->getKey());
      ptrStk.push(tmp[i]);
    }
  }
}

template<typename T>
void Traversal<T>::topDownTraversal(Node<T> *root, CutoffWorker<T> *worker, State *state){
  stack<Node<T>*> ptrStk;
  //CkPrintf("(%d) VISIT node %ld\n", CkMyPe(), getKey());
  int ret = worker->work(root);
  if(ret == 0) return;

  ptrStk.push(root); 
  //CkPrintf("(%d) PUSH node %ld\n", CkMyPe(), getKey());

  Node<T> *tmp[BRANCH_FACTOR];
  int numChildrenToExpand;
  //tmp.reserve(2);

  while(!ptrStk.empty()){
    //pair<Node<T>*,int> pr = ptrStk.top();
    Node<T> *t_parent = ptrStk.top();
    Node<T> *t_children = t_parent->getChildren();
    int t_numChildren = t_parent->getNumChildren();
    //CkPrintf("(%d) POP node %ld\n", CkMyPe(), t_parent->getKey());
    ptrStk.pop();

    //tmp.length() = 0;
    numChildrenToExpand = 0;
    if(t_numChildren == 0){
      processLeaf(t_parent,worker,state);
      continue;
    }

    for(int i = 0; i < t_numChildren; i++){
      Node<T> *t_node = t_children+i;
      //CkPrintf("(%d) VISIT node %ld\n", CkMyPe(), t_node->getKey());
      ret = worker->work(t_node);
      if(ret > 0){
        //tmp.push_back(t_node);
        tmp[numChildrenToExpand] = t_node;
        numChildrenToExpand++;
      }
    }
    for(int i = numChildrenToExpand-1; i >= 0; i--){
      //CkPrintf("(%d) PUSH child %d key %ld\n", CkMyPe(), i, tmp[i]->getKey());
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
    //pair<Node<T>*,int> cpair = cstk.top();
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

