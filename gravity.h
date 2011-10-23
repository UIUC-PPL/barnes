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

/*
 * grav:
 * Calculate the forces on particles pstart to (not including) pend
 * due to a source point 'mass' (possibly another particle or an entire node)
 * located at 'position'.
 */
inline 
void grav(Particle *pstart, Particle *pend, Real mass, const Vector3D<Real> &position){
  Vector3D<Real> dr;
  Real drsq;
  Real drabs;
  Real phii;
  Real mor3;

  for(Particle *p = pstart; p != pend; p++){
    if(position == p->position) continue;
    dr = position - p->position;
    drsq = dr.lengthSquared();
    drsq += globalParams.epssq;
    drabs = sqrt((double) drsq);
    phii = mass/drabs;
    p->potential -= phii;
    mor3 = phii/drsq;
    p->acceleration += mor3*dr;
  }
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
  grav(particles,particles+numParticles,node->data.moments.totalMass,node->data.moments.cm);
  return req->getNumParticles();
}

/*
 * partBucketForce:
 * Calculate force on bucket 'req' due to source particle 'part'
 */
inline int partBucketForce(ExternalParticle *part, 
			   Node<ForceData> *req){ 

  Particle *particles = req->getParticles();
  int numParticles = req->getNumParticles();
  grav(particles,particles+numParticles,part->mass,part->position);
  return numParticles;
}

#endif
