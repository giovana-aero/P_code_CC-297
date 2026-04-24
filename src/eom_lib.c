#include<math.h>
#include<stdlib.h>
#include"../include/eom_lib.h"

#define pi 3.1415926535897932384626433

/*
- gigiaero, 24/04/2026, 1352 hours
*/
void cosspace(double *x,double xi,double xf,int n){
  double midp = (xf - xi)/2.;
  double th_i = pi/(((double) n) - 1.);
  double th = th_i;

  x[0] = xi;
  for(int i=1;i<n;i++){
      x[i] = xi + midp*(1. - cos(th));
      th += th_i;
  }
}

/*
also makes circles when rx = ry
- gigiaero, 24/04/2026,1635 hours
*/
void ellipse(double *x,double *y,double rx,double ry,int n){
  double *th = malloc(sizeof(double)*n);

  linspace(th,0.,2*pi,n);

  for(int i=0;i<n;i++){
    x[i] = rx*cos(th[i]);
    y[i] = ry*sin(th[i]);
  }
    
  free(th);
}

/*
- gigiaero, 24/04/2026, 1313 hours
*/
void linspace(double *x,double xi,double xf,int n){

x[0] = xi;
double step = (xf - xi)/((double) n - 1);

for(int i=1;i<n;i++)
  x[i] = x[i-1] + step;
}