class Particle : public ExternalParticle {
  Key key;
  Vector3D<Real> velocity = new Vector3D<Real>(0.0);
  Vector3D<Real> acceleration = new Vector3D<Real>(0.0);
  Real potential = 0.0;

  int order;

  // FIXME - how to do this
  boolean operator<=(const Particle &other) const { return key <= other.key; }
  boolean operator>=(const Particle &other) const { return key >= other.key; }
  boolean operator>=(const Key &k) const { return key >= k; }

  void pup(PUPer &p){
    super.pup(p);
    p|velocity;
    p|key;
    p|potential;
  }
};

