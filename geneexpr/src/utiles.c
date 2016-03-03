#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Rrun
#ifdef Rrun
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif

double kdtSTEPdotC (double x[], double y[], int p, double d);
SEXP kdtSTEPdotC2R (SEXP Rx, SEXP Ry, SEXP Rp, SEXP Rd);
double kdtQUADRATICdotC (double x[], double y[], int p, double a);
double quadfunc (double v[], int i, int j, double a);
SEXP kdtQUADRATICdotC2R (SEXP Rx, SEXP Ry, SEXP Rp, SEXP Ra);
double kdtERRORdotC (double x[], double y[], int p, double sigma);
SEXP kdtERRORdotC2R (SEXP Rx, SEXP Ry, SEXP Rp, SEXP Rsigma);

double kdtSTEPdotC (double x[], double y[], int p, double d)
{
    /* suppose [-d,d] treated as equality, where d has to be non-negative */
    
    int i, j, mx = 0, my = 0;
    double s = 0;
  
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < i; j++)
        {
            if ((x[j] - x[i] > d && y[j] - y[i] > d) || (x[i] - x[j] > d && y[i] - y[j] > d)) {
                s = s + 1;
            } else if ((x[j] - x[i] > d && y[i] - y[j] > d) || (x[i] - x[j] > d && y[j] - y[i] > d)) {
                s = s - 1;
            }
            
            if (x[i] - x[j] > d || x[j] - x[i] > d) {
                mx = mx + 1;
            }
            
            if (y[i] - y[j] > d || y[j] - y[i] > d) {
                my = my + 1;
            }
        }
    }
    
    if (mx == 0 || my == 0)
        return 0;
    else
        return (s/sqrt(mx)/sqrt(my));
}

SEXP kdtSTEPdotC2R (SEXP Rx, SEXP Ry, SEXP Rp, SEXP Rd)
{
    Rx = coerceVector(Rx,REALSXP);
    Ry = coerceVector(Ry,REALSXP);
    Rp = coerceVector(Rp,INTSXP);
    Rd = coerceVector(Rd,REALSXP);
    
    SEXP Rres;
    
    PROTECT(Rres = allocVector(REALSXP, 1));
    REAL(Rres)[0] = kdtSTEPdotC (REAL(Rx), REAL(Ry), INTEGER(Rp)[0], REAL(Rd)[0]);
    UNPROTECT(1);
    
    return Rres;

}

double kdtQUADRATICdotC (double x[], double y[], int p, double a)
{
    /* suppose uniform noise ~ [-a,a], where a has to be non-negative */
    
    int i, j;
    double s = 0;
    
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < i; j++)
        {
            s = s + quadfunc(x, i, j, a) * quadfunc(y, i, j, a);
        }
    }
    
    int n0 = p*(p-1)/2;
    
    return (s/n0);
}

double quadfunc (double v[], int i, int j, double a)
{
    double t = (v[i] - v[j])/(2 * a);
    
    if (t >= 1)
        return (+1);
    else if (t <= -1)
        return (-1);
    else if (t >= 0 && t <= 1)
        return (2 * t - t * t);
    else
        return (2 * t + t * t);
}

SEXP kdtQUADRATICdotC2R (SEXP Rx, SEXP Ry, SEXP Rp, SEXP Ra)
{
    Rx = coerceVector(Rx,REALSXP);
    Ry = coerceVector(Ry,REALSXP);
    Rp = coerceVector(Rp,INTSXP);
    Ra = coerceVector(Ra,REALSXP);
    
    SEXP Rres;
    
    PROTECT(Rres = allocVector(REALSXP, 1));
    REAL(Rres)[0] = kdtQUADRATICdotC (REAL(Rx), REAL(Ry), INTEGER(Rp)[0], REAL(Ra)[0]);
    UNPROTECT(1);
    
    return Rres;
    
}


double kdtERRORdotC (double x[], double y[], int p, double sigma)
{
    /* suppose Gaussian noise ~ N(0,sigma^2), where sigma has to be non-negative */
    
    int i, j;
    double s = 0;
    
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < i; j++)
        {
            s = s + erf( (x[i] - x[j])/sigma/2 ) * erf( (y[i] - y[j])/sigma/2 );
        }
    }
    
    int n0 = p*(p-1)/2;
    
    return (s/n0);
}

SEXP kdtERRORdotC2R (SEXP Rx, SEXP Ry, SEXP Rp, SEXP Rsigma)
{
    Rx = coerceVector(Rx,REALSXP);
    Ry = coerceVector(Ry,REALSXP);
    Rp = coerceVector(Rp,INTSXP);
    Rsigma = coerceVector(Rsigma,REALSXP);
    
    SEXP Rres;
    
    PROTECT(Rres = allocVector(REALSXP, 1));
    REAL(Rres)[0] = kdtERRORdotC (REAL(Rx), REAL(Ry), INTEGER(Rp)[0], REAL(Rsigma)[0]);
    UNPROTECT(1);
    
    return Rres;
    
}

