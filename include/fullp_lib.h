#include"./num_methods.h"

#ifndef _lib_f_potential_
#define _lib_f_potential_

void calc_A_metrics(int m,int n,double A1[m][n],double A2[m][n],double A3[m][n],
                    double x[m][n],double y[m][n],double J[m][n]);
double calc_contraU(double dphi_dksi,double dphi_deta,double A1_val,
                    double A2_val);
double calc_contraV(double dphi_dksi,double dphi_deta,double A2_val,
                    double A3_val);
void calc_J(int m,int n,double J[m][n],double x[m][n],double y[m][n]);
double calc_rho(double contraU,double contraV,double dphi_dksi,
                double dphi_deta);
void evaluate_delta_form_fullp(int m,int n,sim_prmtrs *config,char *fname_msh_x,
                               char *fname_msh_y);
double max2(double n1,double n2);
// double mean2(double n1,double n2);
double mean4_j(int m,int n,double rho[m][n],int i,int j);
double nu_switch(int m,int n,double rho[m][n],double C,double contraUV,
                 int i,int j,int axis);
double rho_coeffs(int m,int n,double rho[m][n],double nu,int rs,int i,int j,
                  int axis);
int rs_idx(double contraUV);
void solve_af2_2d_rectangular_fullp(int m,int n,double phi[m][n],double x[m][n],
                                    double y[m][n],sim_prmtrs *config);
#endif