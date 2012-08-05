/*--------------------------------------------------------------------
  
  NAS Parallel Benchmarks 2.3 OpenMP C versions - CG

  This benchmark is an OpenMP C version of the NPB CG code.
  
  The OpenMP C versions are developed by RWCP and derived from the serial
  Fortran versions in "NPB 2.3-serial" developed by NAS.

  Permission to use, copy, distribute and modify this software for any
  purpose with or without fee is hereby granted.
  This software is provided "as is" without express or implied warranty.
  
  Send comments on the OpenMP C versions to pdp-openmp@rwcp.or.jp

  Information on OpenMP activities at RWCP is available at:

           http://pdplab.trc.rwcp.or.jp/pdperf/Omni/
  
  Information on NAS Parallel Benchmarks 2.3 is available at:
  
           http://www.nas.nasa.gov/NAS/NPB/

--------------------------------------------------------------------*/
/*--------------------------------------------------------------------

  Authors: M. Yarrow
           C. Kuszmaul

  OpenMP C version: S. Satoh
  
--------------------------------------------------------------------*/

/*
c---------------------------------------------------------------------
c  Note: please observe that in the routine conj_grad three 
c  implementations of the sparse matrix-vector multiply have
c  been supplied.  The default matrix-vector multiply is not
c  loop unrolled.  The alternate implementations are unrolled
c  to a depth of 2 and unrolled to a depth of 8.  Please
c  experiment with these to find the fastest for your particular
c  architecture.  If reporting timing results, any of these three may
c  be used without penalty.
c---------------------------------------------------------------------
*/

#include <dispatch/dispatch.h>
#include <Accelerate/Accelerate.h>
#include "npb-C.h"
#include "gen-matrix.h"
#define STR(s) XSTR(s)
#define XSTR(s) #s

#define NZ  NA*(NONZER+1)*(NONZER+1)+NA*(NONZER+2)
#define STRIDE (1024*(CACHE_LINE_SIZE/sizeof(double)))
#define DIVIDE (NA/STRIDE)
#define RESIDUE (NA%STRIDE)

/* global variables */

/* common /partit_size/ */
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;

/* common /main_int_mem/ */
static int colidx[NZ+1];    /* colidx[1:NZ] */
static int rowstr[NA+1+1];  /* rowstr[1:NA+1] */
static int iv[2*NA+1+1];    /* iv[1:2*NA+1] */
static int arow[NZ+1];      /* arow[1:NZ] */
static int acol[NZ+1];      /* acol[1:NZ] */

/* common /main_flt_mem/ */
static double v[NA+1+1];    /* v[1:NA+1] */
static double aelt[NZ+1];   /* aelt[1:NZ] */
static double a[NZ+1];      /* a[1:NZ] */
static double x[NA+2+1];    /* x[1:NA+2] */
static double z[NA+2+1];    /* z[1:NA+2] */
static double p[NA+2+1];    /* p[1:NA+2] */
static double q[NA+2+1];    /* q[1:NA+2] */
static double r[NA+2+1];    /* r[1:NA+2] */
static double w[NA+2+1];    /* w[1:NA+2] */

/* common /urando/ */
static double amult;
static double tran;

/* function declarations */
static void conj_grad (int colidx[], int rowstr[], double x[], double z[],
               double a[], double p[], double q[], double r[],
               double w[], double *rnorm);

/*--------------------------------------------------------------------
      program cg
--------------------------------------------------------------------*/

int main(int argc, char **argv) {
    int i, j, k, it;
    int nthreads = 1;
    double zeta;
    double rnorm;
    double norm_temp11;
    double norm_temp12;
    double t, mflops;
    char class;
    boolean verified;
    double zeta_verify_value, epsilon;

    firstrow = 1;
    lastrow  = NA;
    firstcol = 1;
    lastcol  = NA;

    if (NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10.0) {
    class = 'S';
    zeta_verify_value = 8.5971775078648;
    } else if (NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12.0) {
    class = 'W';
    zeta_verify_value = 10.362595087124;
    } else if (NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20.0) {
    class = 'A';
    zeta_verify_value = 17.130235054029;
    } else if (NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60.0) {
    class = 'B';
    zeta_verify_value = 22.712745482631;
    } else if (NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110.0) {
    class = 'C';
    zeta_verify_value = 28.973605592845;
    } else {
    class = 'U';
    }

    printf("\n\n NAS Parallel Benchmarks 2.3 OpenMP C version"
       " - CG Benchmark\n");
    printf(" Size: %10d\n", NA);
    printf(" Iterations: %5d\n", NITER);

    naa = NA;
    nzz = NZ;

/*--------------------------------------------------------------------
c  Initialize random number generator
c-------------------------------------------------------------------*/
    tran    = 314159265.0;
    amult   = 1220703125.0;
    zeta    = randlc( &tran, amult );

/*--------------------------------------------------------------------
c  
c-------------------------------------------------------------------*/
    makea(naa, nzz, a, colidx, rowstr, NONZER,
      firstrow, lastrow, firstcol, lastcol, 
      RCOND, arow, acol, aelt, v, iv, SHIFT);
    
/*---------------------------------------------------------------------
c  Note: as a result of the above call to makea:
c        values of j used in indexing rowstr go from 1 --> lastrow-firstrow+1
c        values of colidx which are col indexes go from firstcol --> lastcol
c        So:
c        Shift the col index vals from actual (firstcol --> lastcol ) 
c        to local, i.e., (1 --> lastcol-firstcol+1)
c---------------------------------------------------------------------*/
// #pragma omp parallel private(it,i,j,k)
// {   
    // #pragma omp for nowait  // unnessary
    // for (j = 1; j <= lastrow - firstrow + 1; j++) {
    //   for (k = rowstr[j]; k < rowstr[j+1]; k++) {
    //     colidx[k] = colidx[k] - firstcol + 1;
    //   }
    // }

    /*--------------------------------------------------------------------
    c  set starting vector to (1, 1, .... 1)
    c-------------------------------------------------------------------*/
    // Should be optimized by clang
    for (i = 0; i < NA+2; i++) {
      x[i] = 1.0;
    }
    zeta  = 0.0;

    /*-------------------------------------------------------------------
    c
    c  Do one iteration untimed to init all code and data page tables
    c  (then reinit, start timing, to niter its)
    c-------------------------------------------------------------------*/

    for (it = 1; it <= 1; it++) {

      /*--------------------------------------------------------------------
      c  The call to the conjugate gradient routine:
      c-------------------------------------------------------------------*/
      conj_grad (colidx, rowstr, x, z, a, p, q, r, w, &rnorm);

      /*--------------------------------------------------------------------
      c  zeta = shift + 1/(x.z)
      c  So, first: (x.z)
      c  Also, find norm of z
      c  So, first: (z.z)
      c-------------------------------------------------------------------*/
      norm_temp11 = 0.0;
      norm_temp12 = 0.0;

      // TODO: reduction
      for (j = 1; j <= lastcol-firstcol+1; j++) {
              norm_temp11 = norm_temp11 + x[j]*z[j];
              norm_temp12 = norm_temp12 + z[j]*z[j];
      }
      norm_temp12 = 1.0 / sqrt( norm_temp12 );

      /*--------------------------------------------------------------------
      c  Normalize z to obtain x
      c-------------------------------------------------------------------*/
      // TODO: parallel
      for (j = 1; j <= lastcol-firstcol+1; j++) {
              x[j] = norm_temp12*z[j];
      }
    
    } /* end of do one iteration untimed */

/*--------------------------------------------------------------------
c  set starting vector to (1, 1, .... 1)
c-------------------------------------------------------------------*/
    for (i = 0; i < NA+2; i++) {
         x[i] = 1.0;
    }
    zeta  = 0.0;

//} /* end parallel */

    timer_clear( 1 );
    timer_start( 1 );

/*--------------------------------------------------------------------
c---->
c  Main Iteration for inverse power method
c---->
c-------------------------------------------------------------------*/

    for (it = 1; it <= NITER; it++) {

      /*--------------------------------------------------------------------
      c  The call to the conjugate gradient routine:
      c-------------------------------------------------------------------*/
      conj_grad(colidx, rowstr, x, z, a, p, q, r, w, &rnorm);

      /*--------------------------------------------------------------------
      c  zeta = shift + 1/(x.z)
      c  So, first: (x.z)
      c  Also, find norm of z
      c  So, first: (z.z)
      c-------------------------------------------------------------------*/

      norm_temp11 = 0.0;
      norm_temp12 = 0.0;

      // TODO: reduction
      for (j = 1; j <= lastcol-firstcol+1; j++) {
        norm_temp11 = norm_temp11 + x[j]*z[j];
        norm_temp12 = norm_temp12 + z[j]*z[j];
      }

      norm_temp12 = 1.0 / sqrt( norm_temp12 );

      zeta = SHIFT + 1.0 / norm_temp11;

      if( it == 1 ) {
        printf("   iteration           ||r||                 zeta\n");
      }
      printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);

      /*--------------------------------------------------------------------
      c  Normalize z to obtain x
      c-------------------------------------------------------------------*/

      // TODO: parallel
      for (j = 1; j <= lastcol-firstcol+1; j++) {
              x[j] = norm_temp12*z[j];
      }
    } /* end of main iter inv pow meth */

#if defined(_OPENMP)
#pragma omp master
    nthreads = omp_get_num_threads();
#endif /* _OPENMP */

    timer_stop( 1 );

/*--------------------------------------------------------------------
c  End of timed section
c-------------------------------------------------------------------*/

    t = timer_read( 1 );

    printf(" Benchmark completed\n");

    epsilon = 1.0e-10;
    if (class != 'U') {
    if (fabs(zeta - zeta_verify_value) <= epsilon) {
            verified = TRUE;
        printf(" VERIFICATION SUCCESSFUL\n");
        printf(" Zeta is    %20.12e\n", zeta);
        printf(" Error is   %20.12e\n", zeta - zeta_verify_value);
    } else {
            verified = FALSE;
        printf(" VERIFICATION FAILED\n");
        printf(" Zeta                %20.12e\n", zeta);
        printf(" The correct zeta is %20.12e\n", zeta_verify_value);
    }
    } else {
    verified = FALSE;
    printf(" Problem size unknown\n");
    printf(" NO VERIFICATION PERFORMED\n");
    }

    if ( t != 0.0 ) {
    mflops = (2.0*NITER*NA)
        * (3.0+(NONZER*(NONZER+1)) + 25.0*(5.0+(NONZER*(NONZER+1))) + 3.0 )
        / t / 1000000.0;
    } else {
    mflops = 0.0;
    }

    c_print_results("CG", class, NA, 0, 0, NITER, nthreads, t, 
            mflops, "          floating point", 
            verified, STR(NPBVERSION), STR(COMPILETIME),
            STR(CS1), STR(CS2), STR(CS3), STR(CS4), STR(CS5), STR(CS6), STR(CS7));
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/
static void conj_grad (
    int colidx[],   /* colidx[1:nzz] */
    int rowstr[],   /* rowstr[1:naa+1] */
    double x[],     /* x[*] */
    double z[],     /* z[*] */
    double a[],     /* a[1:nzz] */
    double p[],     /* p[*] */
    double q[],     /* q[*] */
    double r[],     /* r[*] */
    double w[],     /* w[*] */
    double *rnorm )
/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/
    
/*---------------------------------------------------------------------
c  Floaging point arrays here are named as in NPB1 spec discussion of 
c  CG algorithm
c---------------------------------------------------------------------*/
{
    __block double d, sum, rho, rho0, alpha, beta;
    int cgit, cgitmax = 25;

    dispatch_queue_t c_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT,0);
    dispatch_queue_t s_queue = dispatch_queue_create("org.idryman.sQueue",DISPATCH_QUEUE_SERIAL);

/*--------------------------------------------------------------------
c  Initialize the CG algorithm:
c-------------------------------------------------------------------*/
    rho = 0.0;
    for (size_t j = 0; j < NA+2; j++) {
      q[j] = 0.0;
      z[j] = 0.0;
      r[j] = x[j];
      p[j] = r[j];
      w[j] = 0.0;
    }

/*--------------------------------------------------------------------
c  rho = r.r
c  Now, obtain the norm of r: First, sum squares of r elements locally...
c-------------------------------------------------------------------*/
    // TODO: reduction
//    __block double testRho = 0.0;
//    dispatch_apply(NA,s_queue,^(size_t idx){
//      testRho += x[idx+1]*x[idx+1];
//    });
    /*
    dispatch_async(s_queue,^{
      for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
        rho += x[j]*x[j];
      }
    });
    dispatch_apply(DIVIDE, c_queue, ^(size_t idx){
      size_t j = idx * STRIDE+1;
      size_t j_e = j + STRIDE;
      double sum=0.0;
      do {
        sum += x[j]*x[j];
      } while (++j < j_e);
      dispatch_sync(s_queue,^{
        rho += sum;
      });
    });
    */
    rho = cblas_dnrm2(NA,&(x[1]),1);
    /*
    for (size_t j = 1; j <= lastcol-firstcol+1; j++) {
      rho = rho + x[j]*x[j];
    }
    */

/*--------------------------------------------------------------------
c  The conj grad iteration loop
c-------------------------------------------------------------------*/
    for (cgit = 1; cgit <= cgitmax; cgit++) {
      rho0 = rho;
      d = 0.0;
      rho = 0.0;
      
/*--------------------------------------------------------------------
c  q = A.p
c  The partition submatrix-vector multiply: use workspace w
c---------------------------------------------------------------------
C
C  NOTE: this version of the multiply is actually (slightly: maybe %5) 
C        faster on the sp2 on 16 nodes than is the unrolled-by-2 version 
C        below.   On the Cray t3d, the reverse is true, i.e., the 
C        unrolled-by-two version is some 10% faster.  
C        The unrolled-by-8 version below is significantly faster
C        on the Cray t3d - overall speed of code is 1.5 times faster.
*/

/* rolled version */      
      dispatch_apply (NA, c_queue, ^(size_t idx){
        size_t j = idx+1;
        double sum = 0.0;
        for (size_t k = rowstr[j]; k < rowstr[j+1]; k++){
          sum += a[k]*p[colidx[k]];
        }
        w[j] = sum;
      });
//       for (size_t j = 1; j <= lastrow-firstrow+1; j++) {
//           sum = 0.0;
//           for (size_t k = rowstr[j]; k < rowstr[j+1]; k++) {
//             sum = sum + a[k]*p[colidx[k]];
//           }
//           w[j] = sum;
//       }
/* unrolled-by-two version 
    #pragma omp for private(i,k)
      for (j = 1; j <= lastrow-firstrow+1; j++) {
        int iresidue;
        double sum1, sum2;
        i = rowstr[j]; 
        iresidue = (rowstr[j+1]-i) % 2;
        sum1 = 0.0;
        sum2 = 0.0;
        if (iresidue == 1) sum1 = sum1 + a[i]*p[colidx[i]];
        for (k = i+iresidue; k <= rowstr[j+1]-2; k += 2) {
          sum1 = sum1 + a[k]   * p[colidx[k]];
          sum2 = sum2 + a[k+1] * p[colidx[k+1]];
        }
          w[j] = sum1 + sum2;
      }
*/
/* unrolled-by-8 version
    #pragma omp for private(i,k,sum)
         for (j = 1; j <= lastrow-firstrow+1; j++) {
         int iresidue;
             i = rowstr[j]; 
             iresidue = (rowstr[j+1]-i) % 8;
             sum = 0.0;
             for (k = i; k <= i+iresidue-1; k++) {
                 sum = sum +  a[k] * p[colidx[k]];
             }
             for (k = i+iresidue; k <= rowstr[j+1]-8; k += 8) {
                 sum = sum + a[k  ] * p[colidx[k  ]]
                           + a[k+1] * p[colidx[k+1]]
                           + a[k+2] * p[colidx[k+2]]
                           + a[k+3] * p[colidx[k+3]]
                           + a[k+4] * p[colidx[k+4]]
                           + a[k+5] * p[colidx[k+5]]
                           + a[k+6] * p[colidx[k+6]]
                           + a[k+7] * p[colidx[k+7]];
             }
             w[j] = sum;
         }
    */
        
        
        // Clang
        // Maybe better to use cblas copy?
        for (size_t j = 0; j < lastcol+1; j++) {
          q[j] = w[j];
        }

    /*--------------------------------------------------------------------
    c  Clear w for reuse...
    c-------------------------------------------------------------------*/
        // Clang
        for (size_t j = 0; j < lastcol+1; j++) {
          w[j] = 0.0;
        }

    /*--------------------------------------------------------------------
    c  Obtain p.q
    c-------------------------------------------------------------------*/
        //
        /*
        dispatch_async(s_queue,^{
          for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
            d += p[j]*q[j];
          }
        });
        dispatch_apply(DIVIDE,c_queue,^(size_t idx){
          size_t j = idx * STRIDE + 1;
          size_t j_e = j + STRIDE;
          double sum = 0.0;
          do {
            sum += p[j]*q[j];
          } while (++j < j_e);
          dispatch_sync (s_queue, ^{
            d += sum;
          });
        });
        */
        d = cblas_ddot(NA,&(p[1]),1, &(q[1]),1);
        /*
        for (size_t j = 1; j <= lastcol-firstcol+1; j++) {
                d = d + p[j]*q[j];
        }
        */

    /*--------------------------------------------------------------------
    c  Obtain alpha = rho / (p.q)
    c-------------------------------------------------------------------*/
        alpha = rho0 / d;

    /*--------------------------------------------------------------------
    c  Save a temporary of rho
    c-------------------------------------------------------------------*/
        /*  rho0 = rho;*/

    /*---------------------------------------------------------------------
    c  Obtain z = z + alpha*p
    c  and    r = r - alpha*q
    c---------------------------------------------------------------------*/
        // TODO parallel
        /*
        dispatch_async(s_queue,^{
          for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
            z[j] = z[j] + alpha * p[j];
            r[j] = r[j] - alpha * q[j];
          }
        });
        dispatch_apply(DIVIDE,c_queue,^(size_t idx){
          size_t j = idx * STRIDE + 1;
          size_t j_e = j + STRIDE;
          do {
            z[j] = z[j] + alpha * p[j];
            r[j] = r[j] - alpha * q[j];
          } while (++j < j_e);
        });
        */
        cblas_daxpy(NA,  alpha, &(p[1]), 1,  &(z[1]), 1);
        cblas_daxpy(NA, -alpha, &(q[1]), 1,  &(r[1]), 1);
          
        /*
        for (size_t j = 1; j <= lastcol-firstcol+1; j++) {
                z[j] = z[j] + alpha*p[j];
                r[j] = r[j] - alpha*q[j];
        }
        */
                
    /*---------------------------------------------------------------------
    c  rho = r.r
    c  Now, obtain the norm of r: First, sum squares of r elements locally...
    c---------------------------------------------------------------------*/
        // TODO reduction
        /*
        dispatch_async(s_queue,^{
          for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
            rho += r[j]*r[j];
          }
        });
        dispatch_apply(DIVIDE, c_queue, ^(size_t idx){
          size_t j = idx * STRIDE+1;
          size_t j_e = j + STRIDE;
          double sum=0.0;
          do {
            sum += r[j]*r[j];
          } while (++j < j_e);
          dispatch_sync(s_queue,^{
            rho += sum;
          });
        });
        */
        rho = cblas_dnrm2 (NA, &(r[1]), 1);

        /*
        for (size_t j = 1; j <= lastcol-firstcol+1; j++) {
                rho = rho + r[j]*r[j];
        }
        */

    /*--------------------------------------------------------------------
    c  Obtain beta:
    c-------------------------------------------------------------------*/
        beta = rho / rho0;

    /*--------------------------------------------------------------------
    c  p = r + beta*p
    c-------------------------------------------------------------------*/
        // TODO parallel
        /*
        dispatch_async(s_queue,^{
          for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
            p[j] = r[j] + beta * p[j];
          }
        });
        dispatch_apply(DIVIDE,c_queue,^(size_t idx){
          size_t j = idx * STRIDE + 1;
          size_t j_e = j + STRIDE;
          do { 
            p[j] = r[j] + beta * p[j]; 
          } while(++j < j_e);
        });
        */
        cblas_dscal (NA, beta, &(p[1]), 1);
        cblas_daxpy (NA, 1, &(r[1]), 1,  &(p[1]), 1);
        /*
        for (size_t j = 1; j <= lastcol-firstcol+1; j++) {
                p[j] = r[j] + beta*p[j];
        }
        */
      } /* end of do cgit=1,cgitmax */

    /*---------------------------------------------------------------------
    c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
    c  First, form A.z
    c  The partition submatrix-vector multiply
    c---------------------------------------------------------------------*/
        sum = 0.0;
        
        // TODO: NO unrolled version?
        dispatch_apply(NA,c_queue,^(size_t idx){
          size_t j = idx+1;
          double sum = 0.0;
          for (size_t k = rowstr[j]; k < rowstr[j+1]; k++){
            sum += a[k] * z[colidx[k]];
          }
          w[j] = sum;
        });
        // for (size_t j = 1; j <= lastrow-firstrow+1; j++) {
        //   d = 0.0;
        //   for (size_t k = rowstr[j]; k <= rowstr[j+1]-1; k++) {
        //     d = d + a[k]*z[colidx[k]];
        //   }
        //   w[j] = d;
        // }

        // Clang
        for (size_t j = 0; j < lastcol+1; j++) {
          r[j] = w[j];
        }

    /*--------------------------------------------------------------------
    c  At this point, r contains A.z
    c-------------------------------------------------------------------*/
    // TODO: reduction
        /*
        dispatch_apply(DIVIDE,c_queue,^(size_t idx){
          size_t j = idx * STRIDE + 1;
          size_t j_e = j + STRIDE;
          do {
            w[j] = (x[j] - r[j]);
            j++;
          } while (j<j_e);
          dispatch_async(c_queue,^{
            size_t j = idx * STRIDE + 1;
            double d = 0.0;
            do {
              d += w[j]*w[j];
              j++;
            } while (j<j_e);
            dispatch_async(s_queue,^{
              sum += d;
            });
          });
        });
        dispatch_barrier_sync(s_queue,^{
          for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
            w[j] = (x[j] - r[j]);
          }
        });
        dispatch_barrier_sync(s_queue,^{
          for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
            sum += w[j]*w[j];
          }
        });
        */
        /*
        sum = 0.0;
        dispatch_async(s_queue,^{
          for (size_t j = DIVIDE*STRIDE+1; j < NA+1; j++){
            d = x[j] - r[j];
            sum += d*d;
          }
        });
        dispatch_apply(DIVIDE,c_queue,^(size_t idx){
          size_t j = idx * STRIDE + 1;
          size_t j_e = j + STRIDE;
          double dd = 0.0;
          double d = 0.0;
          do {
            dd = x[j] - r[j];
            d += dd*dd;
          } while (++j < j_e);
          dispatch_sync(s_queue,^{
            sum += d;
          });
        });
        */
        cblas_daxpy (NA, -1, &(x[1]), 1,  &(r[1]), 1);
        sum = cblas_dnrm2 (NA, &(r[1]), 1);


        /*
        for (size_t j = 1; j <= lastcol-firstcol+1; j++) {
          d = x[j] - r[j];
          sum = sum + d*d;
        }
        */
        
        (*rnorm) = sqrt(sum);
     dispatch_release (c_queue);
     dispatch_release (s_queue);
    }

