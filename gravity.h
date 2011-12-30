#ifndef __GRAVITY_H__
#define __GRAVITY_H__

/*
 * CharmBH: gravity.h
 * Contains functions for
 * 1. opening criterion calculation
 * 2. interaction between target bucket of particles and source tree node
 * 3. interaction between a target bucket of particles and source bucket
 */

#include "Node.h"
#include "Vector3D.h"
#include "defines.h"
#include "Parameters.h"

extern Parameters globalParams;

/*
 * openCriterionBucket:
 * Is 'node' close enough to the bucket of particles (bucketNode)
 * that its children must be examined separately?
 */
inline bool
openCriterionBucket(Node<ForceData> *node,
                   Node<ForceData> *bucketNode) {
  Vector3D<Real> dr = node->data.moments.cm - bucketNode->data.moments.cm;
  Real drsq = dr.lengthSquared();
  return (globalParams.tolsq*drsq < node->data.moments.rsq);
}

inline bool
openCriterionParticle(Node<ForceData> *node,
                   Particle *target) {
  Vector3D<Real> dr = node->data.moments.cm - target->position;
  Real drsq = dr.lengthSquared();
  return (globalParams.tolsq*drsq < node->data.moments.rsq);
}

/*
 * grav:
 * Calculate the forces on particles pstart to (not including) pend
 * due to a source point 'mass' (possibly another particle or an entire node)
 * located at 'position'.
 */
inline 
void grav(Particle *p, Real mass, const Vector3D<Real> &position){
  Vector3D<Real> dr;
  Real drsq;
  Real drabs;
  Real phii;
  Real mor3;

  dr = position - p->position;
  drsq = dr.lengthSquared();
  drsq += globalParams.epssq;
  drabs = sqrt((double) drsq);
  phii = mass/drabs;
  p->potential -= phii;
  mor3 = phii/drsq;
  p->acceleration += mor3*dr;
}

/*
 * nodeBucketForce:
 * Calculate force on bucket 'req' due to source tree 'node'
 */
inline
int nodeBucketForce(Node<ForceData> *node, 
		    Node<ForceData> *req){
  Particle *particles = req->getParticles();
  int numParticles = req->getNumParticles();
  for(Particle *p = particles; p != particles+numParticles; p++){
    grav(p,node->data.moments.totalMass,node->data.moments.cm);
  }
  return numParticles;
}

inline
int nodeParticleForce(Node<ForceData> *node, 
		      Particle *target){
  grav(target,node->data.moments.totalMass,node->data.moments.cm);
  return 1;
}

/*
 * partBucketForce:
 * Calculate force on bucket 'req' due to source particle 'part'
 */
inline int partBucketForce(ExternalParticle *part, 
			   Node<ForceData> *req){ 
  Particle *particles = req->getParticles();
  int numParticles = req->getNumParticles();
  int computed = 0;

  for(Particle *p = particles; p != particles+numParticles; p++){
    if(p == part) continue;
    grav(p,part->mass,part->position);
    computed++;
  }
  return computed;
}

inline int particleParticleForce(ExternalParticle *part, 
			         Particle *target){
  if(target == part) return 0;
  grav(target,part->mass,part->position);
  return 1;
}

#endif
