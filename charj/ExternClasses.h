struct TreeRootContainer {
  Node<ForceData> *root;
  int pe;

  TreeRootContainer(Node<ForceData> *r, int pe_) : 
    root(r),
    pe(pe_)
  {}
};
