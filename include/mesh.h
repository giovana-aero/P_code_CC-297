#ifndef _lib_mesh_
#define _lib_mesh_

void save_mesh(int m,int n,double *x,double *y,char *casename);
void uniform_rectangular_mesh(int m,int n,double delta_x,double delta_y,
  double *x,double *y);

#endif