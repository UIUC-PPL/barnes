#ifndef __MESSAGES_H__
#define __MESSAGES_H__

/*
 * CharmBH - Messages.h
 * Definitions of message structures used for various tasks.
 */

#include "Particle.h"
#include "MultipoleMoments.h"
#include "Node.h" 

struct ParticleMsg : public CMessage_ParticleMsg {
  Particle *part;
  int numParticles;
};

/*
 * RangeMsg:
 * During decomposition, the master DataManager (DM) communicates the range 
 * of particles held by each tree piece to all other DMs. It uses the RangeMsg
 * for this purpose.
 */
struct RangeMsg : public CMessage_RangeMsg {
  Key *keys;
  int numTreePieces;
};

/*
 * When a DM requests a bucket of particles from another, it receives a 
 * ParticleReplyMsg in return.
 */
struct ParticleReplyMsg : public CMessage_ParticleReplyMsg {
  Key key;
  ExternalParticle *data;
  int np;
};

/*
 * When a DM requests the subtree under a particular node from another DM, 
 * it receives a NodeReplyMsg in return.
 */
struct NodeReplyMsg : public CMessage_NodeReplyMsg {
  Key key;
  Node<ForceData> *data;
  int nn;
};

/* 
 * In order to be responsive to requests for data from remote DMs,
 * a TreePiece must cooperatively relinquish the processor after doing some
 * traversal work. In order to reschedule traversal work, the TreePiece
 * sends a RescheduleMsg to itself.
 */
struct RescheduleMsg : public CMessage_RescheduleMsg {
};
#endif
