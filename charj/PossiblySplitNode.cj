class PossiblySplitNode<T> {
  Node<T> node = null;
  boolean split = false;

  public PossiblySplitNode(Node<T> nd){
    node = nd;
  }

  Node<T> getNode(){
    return node;
  }

  boolean isSplit(){
    return split;
  }
};
