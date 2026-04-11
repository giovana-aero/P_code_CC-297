#include"./num_methods.h"

#ifndef _lib_bi_air_lib_
#define _lib_bi_air_lib_

typedef struct biconvex_airfoil_physics_mesh_config{
  int ILE,ITE;
  int IMAX,JMAX;
  double XSF,YSF;
  double uinf,t;
  int af_type; // 1->biconvex, 2->symmetrical naca4
}bi_air_phys_mesh;

void apply_down_b_c(int m,int n,double phi[m][n],double *x,double *y,
                    bi_air_phys_mesh *b_a_m);
void biconvex_airfoil_mesh(bi_air_phys_mesh *b_a_m,double *x, double*y);
double bi_air_shape_dx(double t,double x_i);
void bi_air_dirichlet_vals_free(double *x,b_conditions_2d *b_c,double uinf);
void bi_air_dirichlet_vals_wall(double *x,b_conditions_2d *b_c,double uinf,
                                double t);
void build_linear_sys_matrix_cols2(int m,int n,double A[m-2][m-2],double *py1,
                                   double *py2);
// void build_linear_sys_matrix_lines2(int m,int n,double A[n-2][n-2],double *px1,
//                                     double *px2);
void evaluate_delta_form_bi_air(int m,int n,double phi[m][n],double *x,
                                double *y,sim_parameters *config,
                                bi_air_phys_mesh *b_a_m);
void get_cp_bi_air(int m,int n,double cp[m][n],double u[m][n],double v[m][n],
                   double *x,double *y,bi_air_phys_mesh *b_a_m);
void get_cp_bi_air_chord(double *cp,int m,int n,double phi[m][n],double u[m][n],
                         double *x,double *y,bi_air_phys_mesh *b_a_m);
void get_u_v_potential(int m,int n,double phi[m][n],double u[m][n],
                       double v[m][n],double Ve[m][n],double *x,double *y);
double naca4_symm_dx(double t,double x_i);
void set_mesh_prmtrs(int mtype,bi_air_phys_mesh *b_a_mesh);
void solve_adi_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                     double *y,sim_parameters *config,
                                     bi_air_phys_mesh *b_a_m);
void solve_g_seidel_2d_rectangular_bi_air(int m,int n,double phi[m][n],
                                          double *x,double *y,
                                          sim_parameters *config,
                                          bi_air_phys_mesh *b_a_m);
void solve_lgs_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                     double *y,sim_parameters *config,
                                     bi_air_phys_mesh *b_a_m);
void solve_p_jacobi_2d_rectangular_bi_air(int m,int n,double phi[m][n],
                                          double *x,double *y,
                                          sim_parameters *config,
                                          bi_air_phys_mesh *b_a_m);
void solve_slor_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                      double *y,sim_parameters *config,
                                      bi_air_phys_mesh *b_a_m);
void solve_sor_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                     double *y,sim_parameters *config,
                                     bi_air_phys_mesh *b_a_m);
#endif