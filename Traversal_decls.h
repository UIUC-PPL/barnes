#ifndef __TRAVERSAL_DECLS_H__
#define __TRAVERSAL_DECLS_H__

#include "Node.h"
#include "Particle.h"
#include "Worker.h"

#include "State.h"

class DataManager;
template<typename T> class CutoffWorker;

template<typename T>
class Traversal {

  // to request particles/nodes from
  DataManager *dm;

  void processLeaf(Node<T> *node, CutoffWorker<T> *worker, State *state);

  public:

  Traversal() : dm(NULL) {}

  void setDataManager(DataManager *mgr){
    dm = mgr;
  }

  // preorder does not have good cache behavior, because
  // all children are stored together but processed
  // separately
  void preorderTraversal(Node<T> *root, CutoffWorker<T> *worker);
  void topDownTraversal(Node<T> *root, CutoffWorker<T> *worker, State *state);
  void topDownTraversal_local(Node<T> *root, CutoffWorker<T> *worker);
  void postorderTraversal(Node<T> *root, CutoffWorker<T> *worker);

  //sph
  void startSphTraversal(Node<T> *root, CutoffWorker<T> *worker, State* st);
  int bottomUpTraversal(Node<T> *root, CutoffWorker<T> *worker, State *state);
  int topDownTraversalSph(Node<T> *root, CutoffWorker<T> *worker, State *state);
  int processSphLeaf(Node<T> *leaf, CutoffWorker<T> *worker, State *state);

};
#endif // __TRAVERSAL_DECLS_H__
