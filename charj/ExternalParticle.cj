class ExternalParticle {
  Vector3D<Real> position = new Vector3D<Real>(0.0);
  Real mass = 0.0;

  void pup(PUPer p){
    p|position;
    p|mass;
  }
};


