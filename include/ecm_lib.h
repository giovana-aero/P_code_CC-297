#include"./num_methods.h"
#include"./eom_lib.h"

#ifndef _lib_ecm_
#define _lib_ecm_

void calc_A_ecm(int m,int n,double A[m][n],double x[m][n],double y[m][n],
                int TE);
void calc_B_ecm(int m,int n,double A[m][n],double x[m][n],double y[m][n],
                int TE);
void calc_C_ecm(int m,int n,double A[m][n],double x[m][n],double y[m][n],
                int TE);
void evaluate_delta_form_ecm(sim_prmtrs *config,msh_prmtrs *msh,
                             control_prmtrs *c_prmtrs,int init_only);
void initialize_mesh_ecm(int m,int n,double x[m][n],double y[m][n],
                         msh_prmtrs *msh);
void half_ellipse(double *x,double *y,double *prmtrs,int n,int invert_th);
void L_phi_ecm(int m,int n,double L_phi_x[m][n],double L_phi_y[m][n],
               double x[m][n],double y[m][n],double A[m][n],
               double B[m][n],double C[m][n],int TE);
void solve_adi_2d_rectangular_ecm(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config,int TE);
double uniform_scheme_der1_o2_central_ecm(int m,int n,double phi[m][n],int i);
double uniform_scheme_der2_o2_central_ecm(int m,int n,double phi[m][n],int i,
                                          int axis);

#endif