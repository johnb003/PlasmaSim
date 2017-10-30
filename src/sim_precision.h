#ifndef SIM_PRECISION_H
#define SIM_PRECISION_H

#include <vector_types.h>  // cuda vectors

#define MAKE_SIMREAL3(a, b, c) make_double3(a, b, c)
#define MAKE_SIMREAL4(a, b, c, d) make_double4(a, b, c, d)
#define RSQRT(x) rsqrt(x)
typedef double SIMREAL;
typedef double3 SIMREAL3;
typedef double4 SIMREAL4;

// #define MAKE_SIMREAL3(a, b, c) make_float3(a, b, c)
// #define MAKE_SIMREAL4(a, b, c, d) make_float4(a, b, c, d)
// #define RSQRT(x) rsqrtf(x)
// typedef float SIMREAL;
// typedef float3 SIMREAL3;
// typedef float4 SIMREAL4;


#endif // SIM_PRECISION_H