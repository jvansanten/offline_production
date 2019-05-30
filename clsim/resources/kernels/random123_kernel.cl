
#ifndef CLSIM_RNG_TYPE
#define CLSIM_RNG_TYPE philox4x32
#endif

#include <Random123/threefry.h>
#include <Random123/philox.h>

#define CLSIM_RNG_SUB_R(x,suffix) x##suffix
#define CLSIM_RNG_SUB(typ,suffix) CLSIM_RNG_SUB_R(typ,suffix)
typedef CLSIM_RNG_SUB(CLSIM_RNG_TYPE,_key_t) clsim_rng_key_t;
typedef CLSIM_RNG_SUB(CLSIM_RNG_TYPE,_ctr_t) clsim_rng_ctr_t;
#undef CLSIM_RNG_SUB
#undef CLSIM_RNG_SUB_R

// Counter-based RNG as described in http://dl.acm.org/citation.cfm?doid=2063405

// prototypes to make some compilers happy
inline float rand_co(clsim_rng_key_t *rnd_k,clsim_rng_ctr_t *rnd_c);
inline float rand_oc(clsim_rng_key_t *rnd_k,clsim_rng_ctr_t *rnd_c);

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 [0,1) 
//////////////////////////////////////////////////////////////////////////////
inline float rand_co(clsim_rng_key_t *rnd_k,clsim_rng_ctr_t *rnd_c)
{
  uint x = CLSIM_RNG_TYPE(*rnd_c,*rnd_k).v[0];
    // if (rnd_c->v[2]==0)
    // dbg_printf("%u: %u\n",rnd_c->v[3],x);
  rnd_c->v[3]++;
  #ifdef USE_NATIVE_MATH
    return native_divide(convert_float_rtz((uint)(x&0xfffffffful)),(float)0x100000000); // OpenCL - native divide
  #else
    return (convert_float_rtz((uint)(x&0xfffffffful))/(float)0x100000000); // OpenCL
  #endif
}

inline float4 rand_co_float4(const clsim_rng_key_t *rnd_k,clsim_rng_ctr_t *rnd_c)
{
  clsim_rng_ctr_t x = CLSIM_RNG_TYPE(*rnd_c,*rnd_k);
  rnd_c->v[3]++;
  return
#ifdef USE_NATIVE_MATH
    native_divide(convert_float4_rtz((*((uint4*)(&x)))&0xffffffffu), (float)0x100000000);
#else
    convert_float4_rtz((*((uint4*)(&x)))&0xffffffffu)/((float)0x100000000);
#endif

}

//////////////////////////////////////////////////////////////////////////////
//   Generates a random number between 0 and 1 (0,1]
//////////////////////////////////////////////////////////////////////////////
inline float rand_oc(clsim_rng_key_t *rnd_k,clsim_rng_ctr_t *rnd_c)
{
  return 1.0f-rand_co(rnd_k,rnd_c);
}

inline float4 rand_oc_float4(clsim_rng_key_t *rnd_k,clsim_rng_ctr_t *rnd_c)
{
  return 1.0f-rand_co_float4(rnd_k,rnd_c);
} 

// typedefs for later use
#define RNG_ARGS clsim_rng_key_t *rnd_key,clsim_rng_ctr_t *rnd_counter
#define RNG_ARGS_TO_CALL rnd_key,rnd_counter
#define RNG_CALL_UNIFORM_CO rand_co(rnd_key,rnd_counter)
#define RNG_CALL_UNIFORM_OC rand_oc(rnd_key,rnd_counter)
