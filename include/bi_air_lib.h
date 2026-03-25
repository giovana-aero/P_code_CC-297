#include"./num_methods.h"

#ifndef _lib_bi_air_lib_
#define _lib_bi_air_lib_

typedef struct biconvex_airfoil_physics_mesh_config{
  int ILE,ITE;
  int IMAX,JMAX;
  double XSF,YSF;
  double uinf,t;
}bi_air_phys_mesh;

void biconvex_airfoil_mesh(bi_air_phys_mesh *b_a_m,double *x, double*y);
double bi_air_shape(double t,double x_i);
void bi_air_dirichlet_vals_free(double *x,b_conditions_2d *b_c,double uinf);
void bi_air_dirichlet_vals_wall(double *x,b_conditions_2d *b_c,double uinf,
                                double t);
void evaluate_delta_form_bi_air(int m,int n,double phi[m][n],double *x,
                                double *y,sim_parameters *config,
                                bi_air_phys_mesh *b_a_m);
void get_u_v_potential(int m,int n,double phi[m][n],double u[m][n],
                       double v[m][n],double Ve[m][n],double *x,double *y);
void solve_p_jacobi_2d_rectangular_bi_air(int m,int n,double phi[m][n],
                                          double *x,double *y,
                                          sim_parameters *config,
                                          bi_air_phys_mesh *b_a_m);
#endif