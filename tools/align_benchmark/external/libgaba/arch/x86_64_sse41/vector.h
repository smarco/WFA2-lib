
/**
 * @file vector.h
 *
 * @brief header for various vector (SIMD) macros
 */
#ifndef _VECTOR_H_INCLUDED
#define _VECTOR_H_INCLUDED

/**
 * @struct v64_mask_s
 *
 * @brief common 64cell-wide mask type
 */
typedef struct v64_mask_s {
	uint16_t m1;
	uint16_t m2;
	uint16_t m3;
	uint16_t m4;
} v64_mask_t;
typedef struct v64_mask_s v64i8_mask_t;

/**
 * @union v64_mask_u
 */
typedef union v64_mask_u {
	v64_mask_t mask;
	uint64_t all;
} v64_masku_t;
typedef union v64_mask_u v64i8_masku_t;

/**
 * @struct v32_mask_s
 *
 * @brief common 32cell-wide mask type
 */
typedef struct v32_mask_s {
	uint16_t m1;
	uint16_t m2;
} v32_mask_t;
typedef struct v32_mask_s v32i8_mask_t;

/**
 * @union v32_mask_u
 */
typedef union v32_mask_u {
	v32_mask_t mask;
	uint32_t all;
} v32_masku_t;
typedef union v32_mask_u v32i8_masku_t;

/**
 * @struct v16_mask_s
 *
 * @brief common 16cell-wide mask type
 */
typedef struct v16_mask_s {
	uint16_t m1;
} v16_mask_t;
typedef struct v16_mask_s v16i8_mask_t;

/**
 * @union v16_mask_u
 */
typedef union v16_mask_u {
	v16_mask_t mask;
	uint16_t all;
} v16_masku_t;
typedef union v16_mask_u v16i8_masku_t;

/**
 * abstract vector types
 *
 * v2i32_t, v2i64_t for pair of 32-bit, 64-bit signed integers. Mainly for
 * a pair of coordinates. Conversion between the two types are provided.
 *
 * v16i8_t is a unit vector for substitution matrices and gap vectors.
 * Broadcast to v16i8_t and v32i8_t are provided.
 *
 * v32i8_t is a unit vector for small differences in banded alignment. v16i8_t
 * vector can be broadcasted to high and low 16 elements of v32i8_t. It can
 * also expanded to v32i16_t.
 *
 * v32i16_t is for middle differences in banded alignment. It can be converted
 * from v32i8_t
 */
#include "v2i32.h"
#include "v2i64.h"
#include "v16i8.h"
#include "v16i16.h"
#include "v32i8.h"
#include "v32i16.h"
#include "v64i8.h"
#include "v64i16.h"

/* conversion and cast between vector types */
#define _from_v16i8_v64i8(x)	(v64i8_t){ (x).v1, (x).v1, (x).v1, (x).v1 }
#define _from_v32i8_v64i8(x)	(v64i8_t){ (x).v1, (x).v2, (x).v1, (x).v2 }
#define _from_v64i8_v64i8(x)	(v64i8_t){ (x).v1, (x).v2, (x).v3, (x).v4 }

#define _from_v16i8_v32i8(x)	(v32i8_t){ (x).v1, (x).v1 }
#define _from_v32i8_v32i8(x)	(v32i8_t){ (x).v1, (x).v2 }
#define _from_v64i8_v32i8(x)	(v32i8_t){ (x).v1, (x).v2 }

#define _from_v16i8_v16i8(x)	(v16i8_t){ (x).v1 }
#define _from_v32i8_v16i8(x)	(v16i8_t){ (x).v1 }
#define _from_v64i8_v16i8(x)	(v16i8_t){ (x).v1 }

/* reversed alias */
#define _to_v64i8_v16i8(x)		(v32i8_t){ (x).v1, (x).v1, (x).v1, (x).v1 }
#define _to_v64i8_v32i8(x)		(v32i8_t){ (x).v1, (x).v2, (x).v1, (x).v2 }
#define _to_v64i8_v64i8(x)		(v32i8_t){ (x).v1, (x).v2, (x).v3, (x).v4 }

#define _to_v32i8_v16i8(x)		(v32i8_t){ (x).v1, (x).v1 }
#define _to_v32i8_v32i8(x)		(v32i8_t){ (x).v1, (x).v2 }
#define _to_v32i8_v64i8(x)		(v32i8_t){ (x).v1, (x).v2 }

#define _to_v16i8_v16i8(x)		(v16i8_t){ (x).v1 }
#define _to_v16i8_v32i8(x)		(v16i8_t){ (x).v1 }
#define _to_v16i8_v64i8(x)		(v16i8_t){ (x).v1 }

#define _cast_v2i64_v2i32(x)	(v2i32_t){ (x).v1 }
#define _cast_v2i32_v2i64(x)	(v2i64_t){ (x).v1 }

#endif /* _VECTOR_H_INCLUDED */
/**
 * end of vector.h
 */
