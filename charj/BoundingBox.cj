/*
 * BoundingBox:
 * Used to calculate the bounding box of all particles in the
 * simulation universe. It also keeps track of particle energy,
 * to see whether there is a drift in the total system energy 
 * through the simulation.
 */
class BoundingBox {
  OrientedBox<double> box;
  int numParticles;
  Real pe;
  Real ke;
  Real mass;

  public BoundingBox(){
    reset();
  }

  public void reset(){
    numParticles = 0;
    box.reset();
    pe = 0.0;
    ke = 0.0;
    mass = 0.0;
  }

  public void grow(Vector3D<Real> v){
    box.grow(v);
  }

  /*
   * This method is called when performing a reduction over 
   * BoundingBox's. It subsumes the bounding box of the 'other'
   * and accumulates its energy in its own. If a PE has no
   * particles, its contributions are not counted.
   */
  public void grow(BoundingBox other){
    if(other.numParticles == 0) return;
    if(numParticles == 0){
      this = other;
    }
    else{
      box.grow(other.box);
      numParticles += other.numParticles;
      pe += other.pe;
      ke += other.ke;
      mass += other.mass;
    }
  }

  void expand(Real pad){
    box.greater_corner = box.greater_corner*pad+box.greater_corner;
    box.lesser_corner = box.lesser_corner-pad*box.lesser_corner;
  }

  void pup(PUPer p){
    p|box;
    p|numParticles;
    p|pe;
    p|ke;
    p|mass;
  }

};


