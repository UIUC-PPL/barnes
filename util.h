#ifndef __UTIL_H__
#define __UTIL_H__

struct ForceData;
template<typename T>
class Node;

int log2ceil(int n);
Key mssb64(Key x);
int mssb64_pos(Key x);

void findSplitters(Particle *particles, int start, int end, int *splitters, Key childKey, int childDepth);

typedef bool (*BinarySearchGEFn)(void *, void *);

template<typename KEY_TYPE, typename OBJ_TYPE>
int binary_search_ge(const KEY_TYPE &check, const OBJ_TYPE *particles, int start, int end){
  int lo = start;
  int hi = end;
  int mid;
  while(lo < hi){
    mid = lo+((hi-lo)>>1);
    if(particles[mid] >= check){
      hi = mid;
    }
    else{
      lo = mid+1;
    }
  }
  return lo;
}

//template<typename T>
void getMomentsFromParticles(Node<ForceData> *bucket);
//template<typename T>
void getMomentsFromChildren(Node<ForceData> *node);
void getBoundingBoxFromParticles(Node<ForceData> *node);
void getBoundingBoxFromChildren(Node<ForceData> *node);
#endif
