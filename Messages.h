#ifndef __MESSAGES_H__
#define __MESSAGES_H__

#include "Particle.h"
#include "MultipoleMoments.h"
#include "Node.h" 

struct RangeMsg : public CMessage_RangeMsg {
  Key *keys;
  int numTreePieces;
};

struct ParticleReplyMsg : public CMessage_ParticleReplyMsg {
  Key key;
  ExternalParticle *data;
  int np;
};

struct NodeReplyMsg : public CMessage_NodeReplyMsg {
  Key key;
  Node<ForceData> *data;
  int nn;
};

struct RescheduleMsg : public CMessage_RescheduleMsg {
};
#endif
