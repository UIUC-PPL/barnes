#include "defines.h"
#include "Particle.h"
#include "util.h"
#include "Node.h"
#include "MultipoleMoments.h"
#include "TreePiece.h"

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

