#pragma once

#include <assert.h> // assert


#ifdef __NVCC__ 
  extern "C" {
  #include "jurassic.h" // ...
  }

  #define UNROLL _Pragma("unroll")
  #define __ext_inline__ inline
  #define restrict __restrict
#else

  #include "jurassic.h" // ...
  
  #define UNROLL
  // Workaround if no NVIDIA compiler is available: __host__ and __device__ are simply ignored
  #define __host__
  #define __device__
  #define __global__
  #define __ext_inline__ extern inline
#endif

#define ptr const restrict

	__host__ __ext_inline__
  void get_tbl() {
      printf("I am at jurassic gpu jr_common.h, calling init_tbl()\n");
		   init_tbl(); //ako ostavim ovak se pozove funkcija init_tbl iz jurassic-scattera-a (suprotno od onoga sto zelim)
        /* ispis programa:
         *
         * CPU printing vector... 1 2 3 4 5 6 7 
         * CPU printing vector... 2 4 6 8 10 12 14 
         * CPU multiplication...
         * CPU printing vector... 2 8 18 32 50 72 98 
         * CPU printing vector... 2 4 6 8 10 12 14 
         * GPU addition...
         * I am at add_vector function in sum_vector.cu, calling get_tbl()
         * I am at jurassic gpu jr_common.h, calling init_tbl()
         * jurassic scatter init table!
         * CPU printing vector... 4 12 24 40 60 84 112 
         * CPU printing vector... 2 4 6 8 10 12 14 
         *
        */




      // init_tbl_added(); //ako ovo odkomentiram se ne kompajlira!
          /* ispis kompajlera:
           * p/software/jusuf/stages/2020/software/binutils/2.34-GCCcore-9.3.0/bin/ld.gold:
           * error: ../jurassic-gpu/libutil.a(jurassic.o): multiple definition of 'init_tbl'
           * /p/software/jusuf/stages/2020/software/binutils/2.34-GCCcore-9.3.0/bin/ld.gold: /tmp/ccJsokLb.o: previous definition here
           * collect2: error: ld returned 1 exit status 
          */
	} // get_tbl
