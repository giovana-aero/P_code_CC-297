#include<stdio.h>
#include<stdlib.h>
#include"../include/1d_arrays.h"
#include"../include/eom_lib.h"
#include"../include/mesh.h"

int main(){

  int n = 20;
  double xi = 1.;
  double xf = 10.;

  double *x = malloc(sizeof(double)*n);
  double *y = malloc(sizeof(double)*n);

  // linspace(x,xi,xf,n);
  // cosspace(x,xi,xf,n);

  // print_1d_array(n,x);

  ellipses(x,y,1.,.5,n);

  save_mesh(n,n,x,y,"test");

  free(x);
  free(y);

  return 0;
}