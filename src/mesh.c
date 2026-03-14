#include"../include/mesh.h"

void uniform_rectangular_mesh(int m,int n,double delta_x,double delta_y,
  double *x,double *y){

x[0] = 0.;
y[0] = 0.;

for(int i=1;i<n;i++)
x[i] = x[i-1] + delta_x;

for(int j=1;j<m;j++)
y[j] = y[j-1] + delta_y;

}