#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "common.h"
#include "Vector3D.h"

#include "pup.h"

struct Particle;
struct ExternalParticle {
  Vector3D<Real> position;
  Real mass;
  ExternalParticle &operator=(const Particle &p);

  void pup(PUP::er &p){
    p|position;
    p|mass;
  }
};

struct Particle : public ExternalParticle {
  Key key;
  Vector3D<Real> velocity;
  Vector3D<Real> acceleration;
  Real potential;

  Particle() : potential(0.0), acceleration(0.0) { }

  bool operator<=(const Particle &other) const { return key <= other.key; }
  bool operator>=(const Particle &other) const { return key >= other.key; }
  bool operator>=(const Key &k) const {return key >= k; }

  void pup(PUP::er &p){
    ExternalParticle::pup(p);
    p|velocity;
    p|key;
    p|potential;
  }
};


#endif
