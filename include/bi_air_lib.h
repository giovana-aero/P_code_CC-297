#include"./b_conditions.h"

#ifndef _lib_bi_air_lib_
#define _lib_bi_air_lib_

typedef struct biconvex_airfoil_mesh_config{
  int ILE,ITE;
  int IMAX,JMAX;
  double XSF,YSF;
}bi_air_mesh;

/*
- gigiaero, 22/03/2026, 2131 hours
*/
void biconvex_airfoil_mesh(bi_air_mesh *b_a_m,double *x, double*y);
double bi_air_shape(double t,double x_i);
void bi_air_dirichlet_vals_free(double *x,b_conditions_2d *b_c,double uinf);
void bi_air_dirichlet_vals_wall(double *x,b_conditions_2d *b_c,double uinf,
                                double t);

#endif