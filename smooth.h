#ifndef __SMOOTH_H
#define __SMOOTH_H

/*
 * Structure for the priority queue.  
 * And utilities for sph calculations.
 */

#include <queue>
#include "State.h"
#include "Parameters.h"

#include <math.h>
#include <fstream>
using namespace std;

extern Parameters globalParams;

//Object for priority queue entry.
class pqSmoothNode
{
 public:
    double fKey;	// distance^2 -> place in priority queue
    Vector3D<double> dx; // displacement of this particle
    ExternalParticle *p; // pointer to rest of particle data
    
    inline bool operator<(const pqSmoothNode& n) const {
	return fKey < n.fKey;
	}
    };

	
//Object to bookkeep a Bucket Smooth Walk.
class NearNeighborState: public State {
public:
    CkVec<pqSmoothNode> *Qs; //priority queue for each particle
    int nParticlesPending;
    int mynParts; 
    bool started;
    int nSmooth;
    
    NearNeighborState() {
        Qs = new CkVec<pqSmoothNode>[2];
	mynParts = 0; 
        }

    NearNeighborState(int nParts, int nS) {
        Qs = new CkVec<pqSmoothNode>[nParts+2];
	for(int k=0; k<mynParts; k++){
	  Qs[k].reserve(globalParams.nSmooth);
	}
	mynParts = nParts; 
	nParticlesPending = -1;
	started = true;
	nSmooth = nS;
    }
    void reset(int nParts, int nS){
        Qs = new CkVec<pqSmoothNode>[nParts+2];
	for(int k=0; k<mynParts; k++){
	  Qs[k].reserve(globalParams.nSmooth);
	}
	mynParts = nParts; 
	nParticlesPending = -1;
	started = true;
	nSmooth = nS;
    }
    void finishBucketSmooth(int iBucket, TreePiece *tp);
    ~NearNeighborState() {
	delete [] Qs; 
        }

    void startParticle(){
    	nParticlesPending++;
    }
    double getCurrentRadius(){
    	CkVec<pqSmoothNode> &Q = Qs[nParticlesPending];
	return Q[0].fKey;
    }
    int compareAddParticle(ExternalParticle* p_new, Node<ForceData>* node){
    	//CkPrintf("compareAddParticle %d \n..", nParticlesPending);
	Particle* particles = node->getParticles();
    	for(int j = 0; j <= node->getNumParticles(); ++j) {
	  Particle* p_current = particles+j;
	  Vector3D<double> dr = p_current->position - p_new->position;
    	  CkVec<pqSmoothNode> &Q = Qs[j+node->getOwnerStart()];
	  pqSmoothNode pqNew;
	  pqNew.fKey = dr.lengthSquared();
	  pqNew.dx = dr;
	  pqNew.p = p_new;
	  if(Q.size() < nSmooth){
	  	Q.push_back(pqNew);
	  }
	  else{//Q.size() >= nSmooth
	      //pop the largest particle, add the new one
	      double rOld2 = Q[0].fKey; // Ball radius
	      if(rOld2 > dr.lengthSquared()){
	          std::pop_heap(&(Q[0]) + 0, &(Q[0]) + nSmooth);
	          Q.resize(Q.size()-1);       // pop if list is full
	      	Q.push_back(pqNew);
	      	std::push_heap(&(Q[0]) + 0, &(Q[0]) + Q.size());
	      }
	  }
	}
	return node->getNumParticles();
    }
    void compareAddParticle(Particle* p_new, Particle* p_current){
    	//CkPrintf("compareAddParticle %d \n..", nParticlesPending);
	Vector3D<double> dr = p_current->position - p_new->position;
    	CkVec<pqSmoothNode> &Q = Qs[nParticlesPending];
	pqSmoothNode pqNew;
	pqNew.fKey = dr.lengthSquared();
	pqNew.dx = dr;
	pqNew.p = p_new;
	if(Q.size() < nSmooth){
		Q.push_back(pqNew);
	}
	else{//Q.size() >= nSmooth
	    //pop the largest particle, add the new one
	    double rOld2 = Q[0].fKey; // Ball radius
	    if(rOld2 > dr.lengthSquared()){
	        std::pop_heap(&(Q[0]) + 0, &(Q[0]) + nSmooth);
	        Q.resize(Q.size()-1);       // pop if list is full
	    	Q.push_back(pqNew);
	    	std::push_heap(&(Q[0]) + 0, &(Q[0]) + Q.size());
	    }
	}
    }
    void printBallRadius(int iter){
    	CkPrintf("Writing the logs.\n");
    	char logFile[100];
	snprintf(logFile, sizeof(logFile), "radius.log.%d", iter);
    	FILE *f = fopen(logFile, "a");
    	//for each particle, print the ball radius
    	for(int i=0; i<mynParts; i++){
		CkVec<pqSmoothNode> &Q = Qs[i];
		fprintf(f, "%f\n", Q[0].fKey);
		//fprintf(f, "%f %f %f %f\n", Q[0].fKey);
	}
	fclose(f);
    }
    void calculateDensity(){
    	CkPrintf("Calculating the density of %d particles.\n", mynParts);
	double startTime = CmiWallTimer();
    	//for each particle calculate the density by diving the total masses 
	//of the nearest particles to the total volume of the ball
	//just calculate the density and do nothing else
    	for(int i=0; i<mynParts; i++){
		CkVec<pqSmoothNode> &Q = Qs[i];
		double total_mass = 0;
    		for(int j=0; j<globalParams.nSmooth; j++){
    		CkPrintf("%d . %d\n", i, j);
			total_mass += Q[j].p->mass;
		}
		double density = total_mass / (M_PI * Q[0].fKey);
   	}
	double totalTime = CmiWallTimer() - startTime;
    	CkPrintf("Density calculations are finished, took %d.\n", totalTime);
    }
};
/*
 * Return true if a sphere centered at pos with radius^2 of rsq
 * is within the bounding box of node.
 * radius^2 is used for efficiency.
 */
static inline bool
ballWithinBounds(OrientedBox<double>& box, Vector3D<double> pos, double rsq)
{
    double dsq = 0.0;
    double delta;
    
    if((delta = box.lesser_corner.x - pos.x) > 0)
	dsq = delta * delta;
    if(rsq < dsq)
	return false;
    else if((delta = pos.x - box.greater_corner.x) > 0)
	dsq = delta * delta;
    if(rsq < dsq)
	return false;
    if((delta = box.lesser_corner.y - pos.y) > 0)
	dsq = delta * delta;
    if(rsq < dsq)
	return false;
    else if((delta = pos.y - box.greater_corner.y) > 0)
	dsq = delta * delta;
    if(rsq < dsq)
	return false;
    if((delta = box.lesser_corner.z - pos.z) > 0)
	dsq = delta * delta;
    if(rsq < dsq)
	return false;
    else if((delta = pos.z - box.greater_corner.z) > 0)
	dsq = delta * delta;
    return (dsq <= rsq);
    }


/*
 * Return true if a sphere centered at pos with radius^2 of rsq
 * intersects the bounding box of node.
 * radius^2 is used for efficiency.
 */
static inline bool
intersect(OrientedBox<Real> box, Vector3D<Real> pos, float rsq)
{
    float dsq = 0.0;
    float delta;
    
    if((delta = box.lesser_corner.x - pos.x) > 0)
	dsq += delta * delta;
    else if((delta = pos.x - box.greater_corner.x) > 0)
	dsq += delta * delta;
    if(rsq < dsq)
	return false;
    if((delta = box.lesser_corner.y - pos.y) > 0)
	dsq += delta * delta;
    else if((delta = pos.y - box.greater_corner.y) > 0)
	dsq += delta * delta;
    if(rsq < dsq)
	return false;
    if((delta = box.lesser_corner.z - pos.z) > 0)
	dsq += delta * delta;
    else if((delta = pos.z - box.greater_corner.z) > 0)
	dsq += delta * delta;
    return (dsq <= rsq);
    }

/*
inline bool
ballsIntersect(TreePiece *ownerTP,
               Node<ForceData> *node, // Node to test
	       State *state){

    Node<ForceData> *myNode = (Node<ForceData> *) computeEntity;
    Particle *particles = ownerTP->getParticles();
    NearNeighborState *nstate = (NearNeighborState *)state;

    double rBucket = myNode->sizeSm + myNode->fKeyMax;
    if(!intersect(node->boundingBox, myNode->centerSm - offset,rBucket*rBucket))
	return 0;
    for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
	double r2 = nstate->Qs[j][0].fKey; // Ball radius^2
	if(intersect(node->boundingBox, particles[j].position - offset, r2)) {
	    return 1;									
	}
    }										
    return 0;
}
*/
#endif
