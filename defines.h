#ifndef __DEFINES_H__
#define __DEFINES_H__

#include "limits.h"
#include "charm++.h"

typedef CmiUInt8 Key;

#define TREE_KEY_BITS (sizeof(Key)*CHAR_BIT)

#define PREAMBLE_INTS 2 
#define PREAMBLE_REALS 1
#define PREAMBLE_SIZE ((PREAMBLE_INTS)*sizeof(int)+(PREAMBLE_REALS)*sizeof(Real))
#define REALS_PER_PARTICLE 8
#define SIZE_PER_PARTICLE ((REALS_PER_PARTICLE)*sizeof(Real))

typedef float Real;

#define BITS_PER_DIM 21
#define BOXES_PER_DIM (1<<(BITS_PER_DIM))

#define NDIMS 3

#define DECOMP_TOLERANCE 1.2
#define BUCKET_TOLERANCE 1.2

#define DIV2(x) ((x)>>1)
#define EVEN(x) ((DIV2(x)<<1)==(x))

#define BRANCH_FACTOR 2
#define LOG_BRANCH_FACTOR 1

const double opening_geometry_factor = 2 / sqrt(3.0);

#define COSMO_CONST(x) (x)
#define CONVERT_TO_COSMO_TYPE

//#define RRDEBUG CkPrintf
#define RRDEBUG /* empty */

//#define TB_DEBUG CkPrintf
#define TB_DEBUG /* empty */

#endif
