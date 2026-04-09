#include<stdio.h>
#include<stdlib.h>
#include"../../include/1d_arrays.h"
#include"../../include/bi_air_lib.h"

int main(){

  bi_air_phys_mesh b_a_mesh;
  // // ----------------------------------------------- original configuration
  b_a_mesh.ILE = 11;   // Leading edge
  b_a_mesh.ITE = 31;   // Trailing edge
  b_a_mesh.IMAX = 41;  // Number of points along x
  b_a_mesh.JMAX = 12;  // Number of points along y
  b_a_mesh.XSF = 1.25; // Stretching factor, x
  b_a_mesh.YSF = 1.25; // Stretching factor, y

  double t = 0.05;

  double *x = malloc(sizeof(double)*b_a_mesh.IMAX);
  double *y = malloc(sizeof(double)*b_a_mesh.JMAX);

  biconvex_airfoil_mesh(&b_a_mesh,x,y);

  for(int i=b_a_mesh.ILE-1,idx=0;i<b_a_mesh.ITE;i++,idx++)
    printf("x[%d] = %f | %f\n",i,x[i],bi_air_shape(t,x[i]));

  free(x);
  free(y);

  return 0;
}