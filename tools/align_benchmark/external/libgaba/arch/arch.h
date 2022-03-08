
/**
 * @file arch.h
 */
#ifndef _ARCH_H_INCLUDED
#define _ARCH_H_INCLUDED

#include <stdint.h>

/**
 * define capability flags for all processors
 */
#define ARCH_CAP_SSE41			( 0x01 )
#define ARCH_CAP_AVX2			( 0x04 )
#define ARCH_CAP_NEON			( 0x10 )
#define ARCH_CAP_ALTIVEC		( 0x20 )


#ifdef __x86_64__
#  if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
/* the compiler is gcc, not clang nor icc */
#    define _ARCH_GCC_VERSION	( __GNUC__ * 100 + __GNUC_MINOR__ * 10 + __GNUC_PATCHLEVEL__ )
#  endif
/* import architectured depentent stuffs and SIMD vectors */
#  if defined(__AVX2__)
#    include "x86_64_avx2/arch_util.h"
#    include "x86_64_avx2/vector.h"
#  elif defined(__SSE4_1__)
#    include "x86_64_sse41/arch_util.h"
#    include "x86_64_sse41/vector.h"
#  elif !defined(ARCH_CAP)		/* arch.h will be included without SIMD flag to check capability */
#    error "No SIMD instruction set enabled. Check if SSE4.1 or AVX2 instructions are available and add `-msse4.1' or `-mavx2' to CFLAGS."
#  endif

/* map reverse-complement sequence out of the canonical-formed address */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

/**
 * @macro arch_cap
 * @brief returns architecture capability flag
 *
 * CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1 &&
 * CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1 &&
 * CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1
 */
#define arch_cap() ({ \
	uint32_t eax, ebx, ecx, edx; \
	__asm__ volatile("cpuid" \
		: "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx) \
		: "a"(0x01), "c"(0x00)); \
	uint64_t sse4_cap = (ecx & 0x180000) != 0; \
	__asm__ volatile("cpuid" \
		: "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx) \
		: "a"(0x07), "c"(0x00)); \
	uint64_t avx2_cap = (ebx & (0x01<<5)) != 0; \
	uint64_t bmi1_cap = (ebx & (0x01<<3)) != 0; \
	__asm__ volatile("cpuid" \
		: "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx) \
		: "a"(0x80000001), "c"(0x00)); \
	uint64_t lzcnt_cap = (ecx & (0x01<<5)) != 0; \
	((sse4_cap != 0) ? ARCH_CAP_SSE41 : 0) | ((avx2_cap && bmi1_cap && lzcnt_cap) ? ARCH_CAP_AVX2 : 0); \
})

#endif

#ifdef AARCH64

/* use x86_64 default */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

#ifdef PPC64

/* use x86_64 default */
#define GREF_SEQ_LIM			( (uint8_t const *)0x800000000000 )

#endif

/*
 * tests
 */
#if !defined(ARCH_CAP) && (!defined(_ARCH_UTIL_H_INCLUDED) || !defined(_VECTOR_H_INCLUDED))
#  error "No SIMD environment detected. Check CFLAGS."
#endif

#ifndef GREF_SEQ_LIM
#  error "No architecuture detected. Check CFLAGS."
#endif


#endif /* #ifndef _ARCH_H_INCLUDED */
/**
 * end of arch.h
 */
