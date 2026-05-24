#include"./num_methods.h"

#ifndef _lib_eom_
#define _lib_eom_

typedef struct control_parameters{
  int L;
  int M;
  double *al;
  double *bm;
  double *cl;
  double *dm;
  int *ksi_l;
  int *ksi_m;
  int *eta_l;
  int *eta_m;
}control_prmtrs;

typedef struct mesh_parameters{
  int IMAX; // number of points (circunference)
  int JMAX; // number of points (normal)
  double c; // chord length
  /* 1-biconvex, 2-naca4, 3-cst */
  int af_type; 
  /*
  1-[t]
  2-[m,p,t]
  3-[(depends on the order of the bernstein polynomial)]
  */
  double *af_prmtrs;
  double *end_prmtrs; // parameters of the last ellipse [rx,ry,dx,dy]
  /*
  1-uniform (ie, always ellipses)
  2-linspace
  3-cosspace_half
  4-parabolic
  */
  int init_type;
}msh_prmtrs;

void alpha_sequence(double *alpha,int *k,int iter,sim_prmtrs *config);
void alpha_sequence_aH(int m,int n,double x[m][n],double y[m][n],
                       sim_prmtrs *config,int *k);
void calc_A(int m,int n,double A[m][n],double x[m][n],double y[m][n]);
void calc_B(int m,int n,double B[m][n],double x[m][n],double y[m][n]);
void calc_C(int m,int n,double C[m][n],double x[m][n],double y[m][n]);
void calc_D(int m,int n,double D[m][n],double x[m][n],double y[m][n]);
double control_P(control_prmtrs *c,int ksi,int eta);
double control_Q(control_prmtrs *c,int ksi,int eta);
void cosspace(double *x,double xi,double xf,int n,int half);
void cst_airfoil(int n_pts,double *x,double *yu,double *yl,double *prmtrs,
                 double c);
void cst_prmtrs(int af,double *prmtrs);
void ellipse(double *x,double *y,double *prmtrs,int n,int invert_th);
void evaluate_delta_form_eom(sim_prmtrs *config,msh_prmtrs *msh,
                             control_prmtrs *c_prmtrs,int init_only);
void init_af_bi_air(double *x,double *y,double *x_axis,int chord_n,
                    msh_prmtrs *msh);
void init_af_cst(double *x,double *y,double *x_axis,int chord_n,
                 msh_prmtrs *msh);
void init_af_naca4(double *x,double *y,double *x_axis,int chord_n,
                   msh_prmtrs *msh);
void init_type1(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void init_type2(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void init_type3(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void initialize_mesh(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void L_phi_eom(int m,int n,double L_phi_x[m][n],double L_phi_y[m][n],
               double x[m][n],double y[m][n],double A[m][n],
               double B[m][n],double C[m][n],double D[m][n],
               control_prmtrs *c_prmtrs);
void linspace(double *x,double xi,double xf,int n);
void malloc_c_prmtrs(control_prmtrs *c_prmtrs);
double min_physical_spacing(int m,int n,double x[m][n],double y[m][n]);
double max_thickness(int m,int n,double y[m][n]);
void naca4(int n,double *x,double *xu,double *xl,double *yu,double *yl,
           double *prmtrs);
void set_control_prmtrs(int c_type,control_prmtrs *c_prmtrs,msh_prmtrs *msh);
void solve_adi_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config,control_prmtrs *c_prmtrs);
void solve_adi_2d_rectangular_eom_np(int m,int n,double x[m][n],double y[m][n],
                                     sim_prmtrs *config,
                                     control_prmtrs *c_prmtrs);
void solve_af2_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config,control_prmtrs *c_prmtrs);
void solve_slor_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                   sim_prmtrs *config,control_prmtrs *c_prmtrs);
double uniform_scheme_der1_o2_central_prdc_ksi(int m,int n,double phi[m][n],
                                               int j);
double uniform_scheme_der2_o2_central_prdc_ksi(int m,int n,double phi[m][n],
                                               int j,int axis);

#endif