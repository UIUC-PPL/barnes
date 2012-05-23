/*
 * CharmBH: gravity.h
 * Contains functions for
 * 1. opening criterion calculation
 * 2. interaction between target bucket of particles and source tree node
 * 3. interaction between a target bucket of particles and source bucket
 */

class Gravity {
  /*
   * openCriterionBucket:
   * Is 'node' close enough to the bucket of particles (bucketNode)
   * that its children must be examined separately?
   */
  static boolean
    openCriterionBucket(Node<ForceData> node, Node<ForceData> bucketNode) {
      Vector3D<Real> dr = node.data.moments.cm - bucketNode.data.moments.cm;
      Real drsq = dr.lengthSquared();
      return (globalParams.tolsq*drsq < node.data.moments.rsq);
    }

  /*
   * grav:
   * Calculate the forces on particles pstart to (not including) pend
   * due to a source point 'mass' (possibly another particle or an entire node)
   * located at 'position'.
   */
  static void grav(Particle p, Real mass, Vector3D<Real> position, Vector3D<Real> dr){
    Real drsq;
    Real drabs;
    Real phii;
    Real mor3;

    dr = position - p.position;
    drsq = dr.lengthSquared();
    drsq += globalParams.epssq;
    drabs = sqrt((double) drsq);
    phii = mass/drabs;
    p.potential -= phii;
    mor3 = phii/drsq;
    // FIXME += allowed?
    p.acceleration += mor3*dr;
  }

  /*
   * nodeBucketForce:
   * Calculate force on bucket 'targetBucket' due to source tree 'node'
   */
  static int nodeBucketForce(Node<ForceData> sourceNode, Node<ForceData> targetBucket){
    Array<Particle> particles = targetBucket.getParticles();
    Vector3D<Real> dr = new Vector3D<Real>();
    for(int i = 0; i < targetBucket.getNumParticles(); i++){
      grav(particles[i],sourceNode.data.moments.totalMass,sourceNode.data.moments.cm,dr);
    }
    delete dr;
    return targetBucket.getNumParticles();
  }

  /*
   * partBucketForce:
   * Calculate force on bucket 'targetBucket' due to source particle 'part'
   */
  static int partBucketForce(ExternalParticle *sourceParticle, Node<ForceData> *targetBucket){ 
    Array<Particle> targetParticles = targetBucket.getParticles();
    int computed = 0;

    Vector3D<Real> dr = new Vector3D<Real>();
    for(int i = 0; i < particles.size(); i++){
      Particle targetParticle = targetParticles[i];
      if(targetParticle == sourceParticle) continue;
      grav(targetParticle,sourceParticle.mass,sourceParticle.position,dr);
      computed++;
    }
    delete dr;
    return computed;
  }
}
