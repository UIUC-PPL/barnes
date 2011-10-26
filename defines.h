#ifndef __DEFINES_H__
#define __DEFINES_H__

/*
 * CharmBH: defines.h
 * Definitions of a number of constants used in the simulation code.
 */

#include "limits.h"
#include "charm++.h"
#include "common.h"

// Number of bits needed for a node or particle's Key (defaults to 64) 
#define TREE_KEY_BITS (sizeof(Key)*CHAR_BIT)

/*
 * When the floating point coordinates of a particle are integerized, 
 * there are BITS_PER_DIM bits available (out of a total of TREE_KEY_BITS)
 * to hold the integer equivalent of each floating point coordinate.
 */
#define BITS_PER_DIM 21
/*
 * Therefore, each integer coordinate of a particle's position can take one
 * of BOXES_PER_DIM values. 
 */
#define BOXES_PER_DIM (1<<(BITS_PER_DIM))

#define NDIMS 3

/*
 * When doing decomposition or tree building, we allow leaves to have up to
 * *_TOLERANCE times the number of particles specfied by user (through -ppc for 
 * decomposition and -b for tree building)
 */
#define DECOMP_TOLERANCE 1.2
#define BUCKET_TOLERANCE 1.2

/*
 * Quick way to find out whether an integer is even. Avoids the '%' operation.
 */
#define DIV2(x) ((x)>>1)
#define EVEN(x) ((DIV2(x)<<1)==(x))

/*
 * When constructing a tree, each node is either a leaf or has BRANCH_FACTOR children.
 */
#define BRANCH_FACTOR 2
#define LOG_BRANCH_FACTOR 1

#define RRDEBUG /* empty */

/*
 * To assign priorities to different tasks.
 */
#define NUM_PRIORITY_BITS (sizeof(int)*CHAR_BIT)

#define REQUEST_MOMENTS_PRIORITY (-8)
#define RECV_MOMENTS_PRIORITY (-7)
#define REQUEST_NODE_PRIORITY (-6)
#define REQUEST_PARTICLES_PRIORITY (-5)
#define RECV_NODE_PRIORITY (-4)
#define RECV_PARTICLES_PRIORITY (-3)
#define REMOTE_GRAVITY_PRIORITY (-2)
#define LOCAL_GRAVITY_PRIORITY (-1)

#endif
