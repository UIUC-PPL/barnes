class Requestor {
  CutoffWorker<ForceData> worker;
  State state;
  Traversal<ForceData> traversal;
  Node<ForceData> currentBucket;

  Requestor(CutoffWorker<ForceData> w, State s, Traversal<ForceData> t, Node<ForceData> b){
    worker = w;
    state = s;
    traversal = t;
    currentBucket = b;
  }

  Requestor(){
    worker = null;
    state = null;
    traversal = null
    currentBucket = null;
  }
};


