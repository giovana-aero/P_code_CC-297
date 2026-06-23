#include"./num_methods.h"

#ifndef _lib_f_potential_
#define _lib_f_potential_

typedef struct fullp_parameters{
  double alpha; // in degrees
  double Ma;
  double C;
  double beta_sub;
  double beta_super;
  int lift; // Activate formulation to take lift into account
}fullp_prmtrs;

typedef struct fullp_beta_parameters{
  double *L2;
  double *Linf;
  double beta_high;
  double beta_low;
}fp_beta_prmtrs;

void beta_initial_config(sim_prmtrs *config,fp_beta_prmtrs *fpb_prmtrs,
                         double beta);
void beta_switch(double *beta,double beta_sub,double beta_super,int supersonic);
void beta_update(double *beta,double L2_now,double Linf_now,int M,
                 fp_beta_prmtrs *fpb_prmtrs,int iter);
void beta_switch(double *beta,double beta_sub,double beta_super,int supersonic);
double calc_A_metrics(int m,int n,double x[m][n],double y[m][n],int i,int j,
                      int A_num,int op);
double calc_Ai(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
               double C,int i,int j);
double calc_Aj(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
               double C,int i,int j);
double calc_contraUV(double dphi_dksi,double dphi_deta,double A12,double A23);
// double calc_contraV(double dphi_dksi,double dphi_deta,int m,int n,
//                     double A2[m][n],double A3[m][n],int i,int j,int op);
double calc_J(int m,int n,double x[m][n],double y[m][n],int i,int j,int op);
double calc_rho(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
                int i,int j,int op);
void evaluate_delta_form_fullp(int m,int n,sim_prmtrs *config,
                               fullp_prmtrs *fp_prmtrs,char *fname_msh_x,
                               char *fname_msh_y);
double freestream_u(double Ma);
double get_mach(double contraU,double contraV,double dphi_dksi,
                double dphi_deta);
double get_dphi_deta(int m,int n,double phi[m][n],double A2_val,double A3_val,
                     double dphi_dksi,int i,int j,int op);
double get_dphi_dksi(int m,int n,double phi[m][n],int i,int j,int op);
double get_dxy_deta(int m,int n,double xy[m][n],int i,int j,int op);
double get_dxy_dksi(int m,int n,double xy[m][n],int i,int j,int op);
void get_u_v_potential_fullp(int m,int n,double phi[m][n],double x[m][n],
                             double y[m][n],double J[m][n],double u[m][n],
                             double v[m][n],double Ve[m][n]);
void initialize_fullp(int m,int n,double phi[m][n],double x[m][n],
                      double y[m][n],fullp_prmtrs *fp_prmtrs);
void L_phi_fullp(int m,int n,double L_phi[m][n],double phi[m][n],double x[m][n],
                 double y[m][n],double C,int af2_flip);
double L_phi_fullp_der_terms_ih(int m,int n,double phi[m][n],double x[m][n],
                                double y[m][n],double C,int i,int j);
double L_phi_fullp_der_terms_jh(int m,int n,double phi[m][n],double x[m][n],
                                double y[m][n],double C,int i,int j);
// void local_q_dervs(double *dphi_dksi,double *dphi_deta,int m,int n,
//                    double phi[m][n],double A2[m][n],double A3[m][n],
//                    int i,int j,int af2_flip);
double max2(double n1,double n2);
// double mean2(double n1,double n2);
double mean4_j(int m,int n,double rho[m][n-1],int i,int j);
void metric_terms(double *ksix,double *ksiy,double *etax,double *etay,int m,
                  int n,double J[m][n],double x[m][n],double y[m][n],
                  int i,int j);
double norm_L2(int m,int n,double L_phi[m][n]);
double nu_switch(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
                 double C,double contraUV,int i,int j,int axis);
double rho_coeffs(int m,int n,double phi[m][n],double x[m][n],double y[m][n],
                  double nu,int rs,int i,int j,int axis);
int rs_idx(double contraUV);
void three_point_pol2_extrp(int m,int n,double A[m][n],double x[m][n],
                            double y[m][n],int i,int end);
void two_point_linear_extrp(int m,int n,double A[m][n],double x[m][n],
                            double y[m][n],int i,int end);
// void solve_adi_2d_rectangular_fullp(int m,int n,double phi[m][n],double J[m][n],
//                                     double A1[m][n],double A2[m][n],
//                                     double A3[m][n],sim_prmtrs *config,
//                                     fullp_prmtrs *fp_prmtrs);
void solve_af2_2d_rectangular_fullp(int m,int n,double phi[m][n],double x[m][n],double y[m][n],double J[m][n],
                                    double A1[m][n],double A2[m][n],
                                    double A3[m][n],sim_prmtrs *config,
                                    fullp_prmtrs *fp_prmtrs);

#endif