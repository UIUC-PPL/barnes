#ifndef MULTIPOLEMOMENTS_H
#define MULTIPOLEMOMENTS_H

/*
 * CharmBH: MultipoleMoments.h
 * Defines structure for basic (first-order)
 * multipole expansion of a node representing a subvolume
 * of space. It includes the node's total mass, center of mass
 * and distance between the center of mass and the farthest
 * enclosed particle.
 */

#include "pup.h"
#include "Vector3D.h"
#include "defines.h"

class FullTreeNode;
class Particle;

/// A representation of a multipole expansion.
class MultipoleMoments {
public:
	/// A physical size for this multipole expansion, calculated by an external function using some other information
	Real rsq;
	/// The total mass represented by this expansion
	Real totalMass;
	/// The center of mass (zeroth order multipole)
	Vector3D<Real> cm;
        MultipoleMoments() {
          clear();
        }
	
	/// Reset this expansion to nothing
        void clear() {
          rsq = 0;
          totalMass = 0;
          cm.x = cm.y = cm.z = 0;
        }

        void pup(PUP::er &p){
          p | rsq;
          p | totalMass;
          p | cm;
        }
};

#endif //MULTIPOLEMOMENTS_H
