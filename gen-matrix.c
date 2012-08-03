#include "npb-C.h"
#include "gen-matrix.h"
#define NZ  NA*(NONZER+1)*(NONZER+1)+NA*(NONZER+2)

static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;

/* common /urando/ */
static double amult;
static double tran;


/*---------------------------------------------------------------------
c       generate the test problem for benchmark 6
c       makea generates a sparse matrix with a
c       prescribed sparsity distribution
c
c       parameter    type        usage
c
c       input
c
c       n            i           number of cols/rows of matrix
c       nz           i           nonzeros as declared array size
c       rcond        r*8         condition number
c       shift        r*8         main diagonal shift
c
c       output
c
c       a            r*8         array for nonzeros
c       colidx       i           col indices
c       rowstr       i           row pointers
c
c       workspace
c
c       iv, arow, acol i
c       v, aelt        r*8
c---------------------------------------------------------------------*/
static void sparse(double a[], int colidx[], int rowstr[], int n,
           int arow[], int acol[], double aelt[],
           int firstrow, int lastrow,
           double x[], boolean mark[], int nzloc[], int nnza);
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[],
           int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);



void makea(
    int n,
    int nz,
    double a[],     /* a[1:nz] */
    int colidx[],   /* colidx[1:nz] */
    int rowstr[],   /* rowstr[1:n+1] */
    int nonzer,
    int firstrow,
    int lastrow,
    int firstcol,
    int lastcol,
    double rcond,
    int arow[],     /* arow[1:nz] */
    int acol[],     /* acol[1:nz] */
    double aelt[],  /* aelt[1:nz] */
    double v[],     /* v[1:n+1] */
    int iv[],       /* iv[1:2*n+1] */
    double shift )
{
    int i, nnza, iouter, ivelt, ivelt1, irow, nzv;

/*--------------------------------------------------------------------
c      nonzer is approximately  (int(sqrt(nnza /n)));
c-------------------------------------------------------------------*/

    firstrow = 1;
    lastrow  = NA;
    firstcol = 1;
    lastcol  = NA;
    naa = NA;
    nzz = NZ;
/*--------------------------------------------------------------------
c  Initialize random number generator
c-------------------------------------------------------------------*/
    tran    = 314159265.0;
    amult   = 1220703125.0;
    randlc( &tran, amult );

    double size, ratio, scale;
    int jcol;

    size = 1.0;
    ratio = pow(rcond, (1.0 / (double)n));
    nnza = 0;

/*---------------------------------------------------------------------
c  Initialize colidx(n+1 .. 2n) to zero.
c  Used by sprnvc to mark nonzero positions
c---------------------------------------------------------------------*/
#pragma omp parallel for 
    for (i = 1; i <= n; i++) {
    colidx[n+i] = 0;
    }
    for (iouter = 1; iouter <= n; iouter++) {
      nzv = nonzer;
      sprnvc(n, nzv, v, iv, &(colidx[0]), &(colidx[n]));
      vecset(n, v, iv, &nzv, iouter, 0.5);
      for (ivelt = 1; ivelt <= nzv; ivelt++) {
          jcol = iv[ivelt];
          if (jcol >= firstcol && jcol <= lastcol) {
          scale = size * v[ivelt];
          for (ivelt1 = 1; ivelt1 <= nzv; ivelt1++) {
              irow = iv[ivelt1];
              if (irow >= firstrow && irow <= lastrow) {
                nnza = nnza + 1;
                if (nnza > nz) {
                    printf("Space for matrix elements exceeded in"
                       " makea\n");
                    printf("nnza, nzmax = %d, %d\n", nnza, nz);
                    printf("iouter = %d\n", iouter);
                    exit(1);
                }
                acol[nnza] = jcol;
                arow[nnza] = irow;
                aelt[nnza] = v[ivelt1] * scale;
              }
          }
          }
      }
      size = size * ratio;
    }

/*---------------------------------------------------------------------
c       ... add the identity * rcond to the generated matrix to bound
c           the smallest eigenvalue from below by rcond
c---------------------------------------------------------------------*/
    for (i = firstrow; i <= lastrow; i++) {
      if (i >= firstcol && i <= lastcol) {
          iouter = n + i;
          nnza = nnza + 1;
          if (nnza > nz) {
            printf("Space for matrix elements exceeded in makea\n");
            printf("nnza, nzmax = %d, %d\n", nnza, nz);
            printf("iouter = %d\n", iouter);
            exit(1);
          }
          acol[nnza] = i;
          arow[nnza] = i;
          aelt[nnza] = rcond - shift;
      }
    }

/*---------------------------------------------------------------------
c       ... make the sparse matrix from list of elements with duplicates
c           (v and iv are used as  workspace)
c---------------------------------------------------------------------*/
    sparse(a, colidx, rowstr, n, arow, acol, aelt,
       firstrow, lastrow, v, &(iv[0]), &(iv[n]), nnza);
}

/*---------------------------------------------------
c       generate a sparse matrix from a list of
c       [col, row, element] tri
c---------------------------------------------------*/
static void sparse(
    double a[],     /* a[1:*] */
    int colidx[],   /* colidx[1:*] */
    int rowstr[],   /* rowstr[1:*] */
    int n,
    int arow[],     /* arow[1:*] */
    int acol[],     /* acol[1:*] */
    double aelt[],  /* aelt[1:*] */
    int firstrow,
    int lastrow,
    double x[],     /* x[1:n] */
    boolean mark[], /* mark[1:n] */
    int nzloc[],    /* nzloc[1:n] */
    int nnza)
/*---------------------------------------------------------------------
c       rows range from firstrow to lastrow
c       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
c---------------------------------------------------------------------*/
{
    int nrows;
    int i, j, jajp1, nza, k, nzrow;
    double xi;

/*--------------------------------------------------------------------
c    how many rows of result
c-------------------------------------------------------------------*/
    nrows = lastrow - firstrow + 1;

/*--------------------------------------------------------------------
c     ...count the number of triples in each row
c-------------------------------------------------------------------*/
#pragma omp parallel for     
    for (j = 1; j <= n; j++) {
      rowstr[j] = 0;
      mark[j] = FALSE;
    }
    rowstr[n+1] = 0;
    
    for (nza = 1; nza <= nnza; nza++) {
      j = (arow[nza] - firstrow + 1) + 1;
      rowstr[j] = rowstr[j] + 1;
    }

    rowstr[1] = 1;
    for (j = 2; j <= nrows+1; j++) {
      rowstr[j] = rowstr[j] + rowstr[j-1];
    }

/*---------------------------------------------------------------------
c     ... rowstr(j) now is the location of the first nonzero
c           of row j of a
c---------------------------------------------------------------------*/
    
/*--------------------------------------------------------------------
c     ... do a bucket sort of the triples on the row index
c-------------------------------------------------------------------*/
    for (nza = 1; nza <= nnza; nza++) {
      j = arow[nza] - firstrow + 1;
      k = rowstr[j];
      a[k] = aelt[nza];
      colidx[k] = acol[nza];
      rowstr[j] = rowstr[j] + 1;
    }

/*--------------------------------------------------------------------
c       ... rowstr(j) now points to the first element of row j+1
c-------------------------------------------------------------------*/
    for (j = nrows; j >= 1; j--) {
      rowstr[j+1] = rowstr[j];
    }
    rowstr[1] = 1;

/*--------------------------------------------------------------------
c       ... generate the actual output rows by adding elements
c-------------------------------------------------------------------*/
    nza = 0;
#pragma omp parallel for    
    for (i = 1; i <= n; i++) {
      x[i] = 0.0;
      mark[i] = FALSE;
    }

    jajp1 = rowstr[1];
    for (j = 1; j <= nrows; j++) {
      nzrow = 0;
      /*--------------------------------------------------------------------
      c          ...loop over the jth row of a
      c-------------------------------------------------------------------*/
      for (k = jajp1; k < rowstr[j+1]; k++) {
            i = colidx[k];
            x[i] = x[i] + a[k];
            if ( mark[i] == FALSE && x[i] != 0.0) {
        mark[i] = TRUE;
        nzrow = nzrow + 1;
        nzloc[nzrow] = i;
        }
      }

      /*--------------------------------------------------------------------
      c          ... extract the nonzeros of this row
      c-------------------------------------------------------------------*/
      for (k = 1; k <= nzrow; k++) {
        i = nzloc[k];
        mark[i] = FALSE;
        xi = x[i];
        x[i] = 0.0;
        if (xi != 0.0) {
          nza = nza + 1;
          a[nza] = xi;
          colidx[nza] = i;
        }
      }
      jajp1 = rowstr[j+1];
      rowstr[j+1] = nza + rowstr[1];
    } // end for (j = 1; j <= nrows; j++)
}

/*---------------------------------------------------------------------
c       generate a sparse n-vector (v, iv)
c       having nzv nonzeros
c
c       mark(i) is set to 1 if position i is nonzero.
c       mark is all zero on entry and is reset to all zero before exit
c       this corrects a performance bug found by John G. Lewis, caused by
c       reinitialization of mark on every one of the n calls to sprnvc
---------------------------------------------------------------------*/
static void sprnvc(
    int n,
    int nz,
    double v[],     /* v[1:*] */
    int iv[],       /* iv[1:*] */
    int nzloc[],    /* nzloc[1:n] */
    int mark[] )    /* mark[1:n] */
{
    int nn1;
    int nzrow, nzv, ii, i;
    double vecelt, vecloc;

    nzv = 0;
    nzrow = 0;
    nn1 = 1;
    do {
    nn1 = 2 * nn1;
    } while (nn1 < n);

/*--------------------------------------------------------------------
c    nn1 is the smallest power of two not less than n
c-------------------------------------------------------------------*/

    while (nzv < nz) {
      vecelt = randlc(&tran, amult);

      /*--------------------------------------------------------------------
      c   generate an integer between 1 and n in a portable manner
      c-------------------------------------------------------------------*/
      vecloc = randlc(&tran, amult);
      i = icnvrt(vecloc, nn1) + 1;
      if (i > n) continue;

      /*--------------------------------------------------------------------
      c  was this integer generated already?
      c-------------------------------------------------------------------*/
      if (mark[i] == 0) {
          mark[i] = 1;
          nzrow = nzrow + 1;
          nzloc[nzrow] = i;
          nzv = nzv + 1;
          v[nzv] = vecelt;
          iv[nzv] = i;
      }
    }

    for (ii = 1; ii <= nzrow; ii++) {
    i = nzloc[ii];
    mark[i] = 0;
    }
}

/*---------------------------------------------------------------------
* scale a double precision number x in (0,1) by a power of 2 and chop it
*---------------------------------------------------------------------*/
static int icnvrt(double x, int ipwr2) {
    return ((int)(ipwr2 * x));
}

/*--------------------------------------------------------------------
c       set ith element of sparse vector (v, iv) with
c       nzv nonzeros to val
c-------------------------------------------------------------------*/
static void vecset(
    int n,
    double v[], /* v[1:*] */
    int iv[],   /* iv[1:*] */
    int *nzv,
    int i,
    double val)
{
    int k;
    boolean set;

    set = FALSE;
    for (k = 1; k <= *nzv; k++) {
    if (iv[k] == i) {
            v[k] = val;
            set  = TRUE;
    }
    }
    if (set == FALSE) {
    *nzv = *nzv + 1;
    v[*nzv] = val;
    iv[*nzv] = i;
    }
}
