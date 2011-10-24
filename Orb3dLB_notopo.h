#ifndef _ORB3DLB_NOTOPO_H_
#define _ORB3DLB_NOTOPO_H_

#include "CentralLB.h"
#include "Orb3dLB_notopo.decl.h"
#include "TaggedVector3D.h"
#include "OrientedBox.h"

#define XDIM 0
#define YDIM 1
#define ZDIM 2
#define NDIMS 3

#define LEFT_PARTITION 0
#define RIGHT_PARTITION 1
#define INVALID_PARTITION -1

void CreateOrb3dLB_notopo();
BaseLB * AllocateOrb3dLB_notopo();

struct OrbObject {
  int partition;
  int lbindex;
  Vector3D<float> centroid;
  OrbObject() : partition(-1), lbindex(-1) {}
  OrbObject(int tag) : partition(-1), lbindex(tag) {}
};

struct Event {
  int owner;
  float load;
  float position;

  Event(float pos, float ld, int o) : 
    position(pos),
    load(ld),
    owner(o)
  {
  }

  Event() : 
    owner(-1),
    load(0.0),
    position(0.0)
  {
  }

  bool operator<=(const Event &e){
    return position <= e.position;
  }

  bool operator>=(const Event &e){
    return position >= e.position;
  }
};



class Orb3dLB_notopo : public CentralLB {
private:
  CkVec<int> *mapping;
  CkVec<int> *from;

  CkVec<OrbObject> tps;
  CkVec<float> procload;
  CkVec<OrientedBox<float> > procbox;

  // things are stored in here before work
  // is ever called.
  TaggedVector3D *tpCentroids;
  CkReductionMsg *tpmsg;
  int nrecvd;
  bool haveTPCentroids;

  int nextProc;

  CmiBool QueryBalanceNow(int step);
  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);
  void orbPartition(CkVec<Event> *events, OrientedBox<float> &box, int procs);
  int partitionRatioLoad(CkVec<Event> &events, float ratio);

public:
  Orb3dLB_notopo(const CkLBOptions &);
  Orb3dLB_notopo(CkMigrateMessage *m):CentralLB(m) { lbname = "Orb3dLB_notopo"; }
  void work(BaseLB::LDStats* stats);
  void receiveCentroids(CkReductionMsg *msg);
};

#endif /* _Orb3dLB_notopo_H_ */
