#include"./num_methods.h"

#ifndef _lib_f_potential_
#define _lib_f_potential_

typedef struct fullp_parameters{
  double alpha;
  double Ma;
  double C;
  int lift; // Activate formulation to take lift into account
}fullp_prmtrs;

void calc_A_metrics(int m,int n,double A1[m][n],double A2[m][n],double A3[m][n],
                    double x[m][n],double y[m][n],double J[m][n]);
double calc_contraU(double dphi_dksi,double dphi_deta,double A1_val,
                    double A2_val);
double calc_contraV(double dphi_dksi,double dphi_deta,double A2_val,
                    double A3_val);
void calc_J(int m,int n,double J[m][n],double x[m][n],double y[m][n]);
void calc_rho(int m,int n,double rho[m][n-1],double phi[m][n],double A1[m][n],
              double A2[m][n],double A3[m][n]);
void evaluate_delta_form_fullp(int m,int n,sim_prmtrs *config,
                               fullp_prmtrs *fp_prmtrs,char *fname_msh_x,
                               char *fname_msh_y);
double freestream_u(double Ma);
void get_u_v_potential_fullp(int m,int n,double phi[m][n],double u[m][n],
                             double v[m][n],double Ve[m][n]);
void initialize_fullp(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
                      fullp_prmtrs *fp_prmtrs);
double L_phi_fullp_der_terms_ih(int m,int n,double phi[m][n],double J[m][n],
                                double A1[m][n],double A2[m][n],
                                double rho[m][n-1],double C,int i,int j);
double L_phi_fullp_der_terms_jh(int m,int n,double phi[m][n],double J[m][n],
                                double A2[m][n],double A3[m][n],
                                double rho[m][n-1],double C,int i,int j);
void L_phi_fullp(int m,int n,double L_phi[m][n],double phi[m][n],
                 double J[m][n],double A1[m][n],double A2[m][n],double A3[m][n],
                 double rho[m][n-1],double C);
double max2(double n1,double n2);
// double mean2(double n1,double n2);
double mean4_j(int m,int n,double rho[m][n-1],int i,int j);
double nu_switch(int m,int n,double rho[m][n-1],double C,double contraUV,
                 int i,int j,int axis);
double rho_coeffs(int m,int n,double rho[m][n-1],double nu,int rs,int i,int j,
                  int axis);
int rs_idx(double contraUV);
void solve_af2_2d_rectangular_fullp(int m,int n,double phi[m][n],double x[m][n],
                                    double y[m][n],sim_prmtrs *config,
                                    fullp_prmtrs *fp_prmtrs);

#endif