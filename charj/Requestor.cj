class Requestor {
  CutoffWorker worker = null;
  State state = null;
  Traversal traversal = null;
  Node currentBucket = null;

  Requestor(CutoffWorker w, State s, Traversal t, Node b){
    worker = w;
    state = s;
    traversal = t;
    currentBucket = b;
  }
};
