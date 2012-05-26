
// sequential functions that are best kept separate
// from charj

void findOrbLB(){
  LBDatabase *lbdb = LBDatabaseObj();
  numLB = lbdb->getNLoadBalancers();
  BaseLB **lb = lbdb->getLoadBalancers();
  haveOrbLB = false;
  for(int i = 0; i < numLB; i++){
    if(string(lb[i]->lbName()) == "Orb3dLB_notopo"){
      orbLBProxy = lb[i]->getGroupID();
      haveOrbLB = true;
      break;
    }
  }
}

void getNumParticles();
void setParameters(CkArgMsg *m, Parameters &params){
  map<string,string> table;
  params.extractParameters(m->argc, m->argv, table); 

  map<string,string>::iterator it = table.find("in");
  if(it == table.end()){
    CkPrintf("[Main] Must specify input particle file:\n");
    usage();
    CkAbort("bad command line arguments\n");
  }

  /* Input file */
  string filename = params.getsparam("in", table);
  params.nchars = filename.length();
  params.filename = new char[params.nchars+1];
  filename.copy(params.filename,params.nchars);
  params.filename[params.nchars] = '\0';

  /* Opening angle */
  params.theta = params.getrparam("theta", DEFAULT_THETA, table);
  CkPrintf("theta: %f\n", params.theta);

  /* Step time */
  params.dtime = params.getrparam("dtime", DEFAULT_DTIME, table);
  CkPrintf("dtime: %f\n", params.dtime);

  /* Half of step time; used in leapfrog "kick" */
  params.dthf = params.dtime/2.0;

  /* Softening parameter */
  params.epssq = params.getrparam("eps", DEFAULT_EPS, table);
  CkPrintf("eps: %f\n", params.epssq);
  params.epssq = params.epssq*params.epssq;

  params.tolsq = params.theta;
  CkPrintf("tol: %f\n", params.tolsq);
  params.tolsq = params.tolsq*params.tolsq;

  /* Particles per chare (number of particles per tree piece) */
  params.ppc = params.getiparam("ppc", DEFAULT_PPC, table); 
  CkPrintf("ppc: %d\n", params.ppc);

  params.ppb = params.getiparam("b", DEFAULT_PPB, table);
  CkPrintf("bucketSize: %d\n", params.ppb);

  /* Number of steps to run the simulation for */
  params.iterations = params.getiparam("killat", DEFAULT_KILLAT, table);
  CkPrintf("killat: %d\n", params.iterations);

  /* 
     When a remote node is requested, the subtree beneath
     it up to "cacheLineSize" is communicated alongwith the
     node itself. This helps to coarsen the grain of communication.
     */
  params.cacheLineSize = params.getiparam("chunkdepth", DEFAULT_CHUNK_DEPTH, table);
  CkPrintf("chunkdepth: %d\n", params.cacheLineSize);

  /*
     When traversing the global tree for each bucket of particles, 
     do "yield" number of buckets before returning to the scheduler
     to look for outstanding request messages, etc.
     */
  params.yieldPeriod = params.getiparam("yield", DEFAULT_YIELD_PERIOD, table);
  CkPrintf("yieldPeriod: %d\n", params.yieldPeriod);

  /*
     Do load balancing every "balancePeriod" iterations.
     */
  params.balancePeriod = params.getiparam("balancePeriod", DEFAULT_BALANCE_PERIOD, table);
  CkPrintf("balancePeriod: %d\n", params.balancePeriod);

  params.doPrintAccel = false;
  if(table.find("printAccel") != table.end()) params.doPrintAccel = true;
  CkPrintf("doPrintAccel: %d\n", params.doPrintAccel);

  params.doPrintTree = false;
  if(table.find("printTree") != table.end()) params.doPrintTree = true;
  CkPrintf("doPrintTree: %d\n", params.doPrintTree);

  /* Number of steps between decompositions */
  params.decompPeriod = params.getiparam("decompPeriod", DEFAULT_DECOMP_PERIOD, table); 
  CkPrintf("decompPeriod: %d\n", params.decompPeriod);

  /* To prevent sending of very large messages after decomposition */
  params.particleMsgMaxSize = params.getiparam("particleMsgMaxSize", DEFAULT_PARTICLE_MSG_MAX_SIZE, table); 
  CkPrintf("particleMsgMaxSize: %d\n", params.particleMsgMaxSize);

  /* Decide several bits in one iteration when decomposing */
  params.decompLevels = params.getiparam("decompLevels", DEFAULT_DECOMP_BITS, table); 
  CkPrintf("decompLevels: %d\n", params.decompLevels);


  /*
     Use the filename obtained previously to read just the 
     number of particles from the input file.
     */
  getNumParticles();

  /* 
     Find the appropriate number of tree pieces to generate.
     */
  it = table.find("p");
  /* 
     Find the appropriate number of tree pieces to generate.
     */
  if(it == table.end()){
    /*
       If the number of pieces is not specified by the user, 
       calculate it from the number of particles per chare and
       the total number of particles.
       */
    params.numTreePieces = CkNumPes()*8;
  }
  else{
    /*
       If the number of pieces is specified by the user, use that.
       */
    params.numTreePieces = atoi(it.second.c_str());
  }

  CkPrintf("tree pieces: %d\n", params.numTreePieces);
}

void getNumParticles(){
  CkPrintf("[Main] file %s\n", globalParams.filename);

  string filename(globalParams.filename);
  Tipsy::TipsyReader r(filename);
  if(!r.status()) {
    cerr << CkMyPe() << ": Fatal: Couldn't open tipsy file! " << filename << endl;
    CkExit();
    return;
  }

  Tipsy::header tipsyHeader = r.getHeader();
  int nTotalParticles = tipsyHeader.nbodies;
  int nTotalSPH = tipsyHeader.nsph;
  int nTotalDark = tipsyHeader.ndark;
  int nTotalStar = tipsyHeader.nstar;

  CkPrintf("[Main] nsph %d ndark %d nstar %d\n", tipsyHeader.nsph, tipsyHeader.ndark, tipsyHeader.nstar);
}

int getTipsyNBodies(Tipsy::TipsyReader *r){
  return r->getHeader().nbodies;
}

void setTipsyDarkParticle(Particle &p, Tipsy::TipsyReader &r){
  Tipsy::dark_particle dp;
  CkAssert(r.getNextDarkParticle(dp)); 
  p.mass = dp.mass;
  p.position = dp.pos;
  p.velocity = dp.vel;
}

void doPrintTree(string name){
  ostringstream oss;
  CkPrintf("[%d] printing tree\n", CkMyPe());
  oss << name << "." << CkMyPe() << "." << iteration << ".dot";
  ofstream ofs(oss.str().c_str());
  ofs << "digraph " << name << "_" << CkMyPe() << "_" << iteration << " {" << endl;
  if(root != NULL) printTree(root,ofs);
  ofs << "}" << endl;
  ofs.close();
}

void printTree(Node nd, ostream os){
  os << nd.getKey() 
    << "[label=\""<< nd.getKey() 
    << "," << nd.getNumParticles() 
    << "," << nd.getOwnerStart() << ":" << nd.getOwnerEnd() 
    << "\\n" << nd.data.moments.cm
    << "\","
    << "style=\"filled\""
    << "color=\"" << NodeTypeColor[nd.getType()] << "\""
    << "]" << endl;
  //if(nd->getOwnerEnd()-1 == nd->getOwnerStart()) return;
  for(int i = 0; i < nd.getNumChildren(); i++){
    Node *child = nd.getChildren()+i;
    CkAssert(child != NULL);
    os << nd.getKey() << " -> " << child.getKey() << endl;
    printTree(child,os);
  }
}

struct ForceData;
class Node;

int log2ceil(int n);
Key mssb64(Key x);
int mssb64_pos(Key x);

void findSplitters(CkVec<Particle> *particles, int start, int end, int *splitters, Key childKey, int childDepth);

template<typename KEY_TYPE, typename OBJ_TYPE>
int binary_search_ge(const KEY_TYPE &check, const OBJ_TYPE *particles, int start, int end){
  int lo = start;
  int hi = end;
  int mid;
  while(lo < hi){
    mid = lo+((hi-lo)>>1);
    if((*particles)[mid] >= check){
      hi = mid;
    }
    else{
      lo = mid+1;
    }
  }
  return lo;
}

int log2ceil(int n){
  int cur = 1;
  int ret = 0;
  while(cur < n){
    cur <<= 1;
    ret++;
  }
  return ret;
}

int mssb64_pos(Key x) {
  int n;

  if (x == Key(0)) return -1;
  n = 0;
  if (x > Key(0x00000000FFFFFFFF)) {n += 32; x >>= 32;}
  if (x > Key(0x000000000000FFFF)) {n += 16; x >>= 16;}
  if (x > Key(0x00000000000000FF)) {n += 8; x >>= 8;}
  if (x > Key(0x000000000000000F)) {n += 4; x >>= 4;}
  if (x > Key(0x0000000000000003)) {n += 2; x >>= 2;}
  if (x > Key(0x0000000000000001)) {n += 1;}
  return n;
}

Key mssb64(Key x)
{
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  x |= (x >> 32);
  return(x & ~(x >> 1));
}

bool CompareKeys(void *a, void *b){
  Key *ownerKey = (Key *)a;
  Key *k = (Key *)b;
  return (*ownerKey >= *k);
}

bool CompareParticleToKey(void *a, void *b){
  Particle *p = (Particle *)a;
  Key *k = (Key *)b;
  return (p->key >= *k);
}

void findSplitters(Particle *particles, int start, int end, int *splitters, Key childKey, int childDepth){
  int nRankBits = LOG_BRANCH_FACTOR;
  /*
  CkPrintf("findSplitters parentKey %lx nRankBits %d startChildKey %lx parentDepth %d childDepth %d\n", 
      parentKey, 
      nRankBits, childKey,
      parentDepth, childDepth);
  */
  // particles of first child always begin at index 0
  splitters[0] = start;
  // for all other children, set the beginning of each 
  // (and, hence, the end of the previous child)
  for(int i = 1; i < BRANCH_FACTOR; i++){
    childKey++;
    Key testKey = (childKey << (TREE_KEY_BITS-(nRankBits*childDepth+1))); 
    //CkPrintf("(%d) findSplitters [%d - %d] start %lx end %lx test %lx\n", CkMyPe(), start, end-1, particles[start].key, particles[end-1].key, testKey);
    int firstGEIdx = binary_search_ge<Key,Particle>(testKey, particles, start, end); 
    splitters[i] = firstGEIdx;
    start = firstGEIdx;
  }
  // end of the last child is always at index end 
  // (which is one position past the last particle)
  splitters[BRANCH_FACTOR] = end;
}

ExternalParticle &ExternalParticle::operator=(const Particle &p){
  this->position = p.position;
  this->mass = p.mass;
  return *this;
}

