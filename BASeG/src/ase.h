/*
 *  ase.h
 *
 *  Created by Wei Sun on 5/25/2010.
 *  Modified by Vasyl Zhabotynsky on 04/05/2011
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"
#include "lbfgsb1.h"

double negLogH0 (int n, double* para, void* ex, SEXP x1);
double negLogH1 (int n, double* para, void* ex, SEXP x1);

void negGradLog (int n, double* para, double* gr, void* ex, SEXP x1);
double negLog (int n, double* para, void* ex, SEXP x1);

void negGradLogH0 (int n, double* para, double* gr, void* ex, SEXP x1);
void negGradLogH1 (int n, double* para, double* gr, void* ex, SEXP x1);


double h3p(double v, double u, double y1, double y2, double mu1, double mu2, double *mean0, double *sigma, double p, int trace);
double h3preal(double v, double u, double y1, double y2, double mu1, double mu2, double *mean0, double *sigma, double p, int trace);
double h3sk(int k, double v, double u, double y1, double y2, double mu1, double mu2, double *mean0, double *sigma, double p, int trace);
double h3skreal(int k, double v, double u, double y1, double y2, double mu1, double mu2, double *mean0, double *sigma, double p, int trace);
double h3bjk(int j, int k, double v, double u, double y1, double y2, double mu1, double mu2, double x1j, double x2j, double *mean0, double *sigma, double p, int trace);

double h4p(double v, double u,  double n1, double n2, double nA1, double nA2, double pi1, double pi2, double *sigma, double p, int trace);
double h4pik(int k, double v, double u, double n1, double n2, double nA1, double nA2, double pi1, double pi2, double *sigma, double p, int trace);
double h4skreal(int k, double v, double u, double n1, double n2, double nA1, double nA2, double pi1, double pi2, double *sigma, double p, int trace);

void ase (int* dims, double* Y1, double* Y2, double* Z, char** output, 
          double* RP_cut, int* cis_only, int* cis_distance, 
          int* eChr, int* ePos, int* mChr, int* mPos,  
          int* trace, int* succeed);
