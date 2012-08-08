#ifndef GEN_MATRIX
#define GEN_MATRIX
void makea(int n, int nz, double a[], int colidx[], int rowstr[],
          int nonzer, int firstrow, int lastrow, int firstcol,
          int lastcol, double rcond, int arow[], int acol[],
          double aelt[], double v[], int iv[], double shift );
#endif
