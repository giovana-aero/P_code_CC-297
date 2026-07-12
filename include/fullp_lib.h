#include"./num_methods.h"

#ifndef _lib_f_potential_
#define _lib_f_potential_

typedef struct fullp_parameters{
  double alpha; // Angle of attack in degrees
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
double calc_Ai(int m,int n,double Ai[m][n],double rho_til[m][n],
               double A1_ih[m][n],double J_ih[m][n]);
double calc_Aj(int m,int n,double Aj[m][n],double rho_bar[m][n],
               double A3_jh[m][n],double J_jh[m][n]);
void beta_switch(double *beta,double beta_sub,double beta_super,int supersonic);
void calc_J_A_metrics(int m,int n,double J[m][n],double A1[m][n],
                      double A2[m][n],double A3[m][n],double x[m][n],
                      double y[m][n],int op);
void calc_contraUV(int m,int n,double contraUV[m][n],double phi[m][n],
                   double A1[m][n],double A2[m][n],double A3[m][n],int op,
                   int axis);
double calc_rho(int m,int n,double rho[m][n],double phi[m][n],double A2[m][n],
                double A3[m][n],double contraU[m][n],double contraV[m][n],
                int op);
void evaluate_delta_form_fullp(int m,int n,sim_prmtrs *config,
                               fullp_prmtrs *fp_prmtrs,char *fname_msh_x,
                               char *fname_msh_y);
double freestream_u(double Ma);
// double get_mach(double contraU,double contraV,double dphi_dksi,
//                 double dphi_deta);
double get_dphi_deta(int m,int n,double phi[m][n],double A2_val,double A3_val,
                     double dphi_dksi,int i,int j,int op);
double get_dphi_dksi(int m,int n,double phi[m][n],int i,int j,int op);
double get_dxy_deta(int m,int n,double xy[m][n],int i,int j,int op);
double get_dxy_dksi(int m,int n,double xy[m][n],int i,int j,int op);
void get_half_meshes(int m,int n,double x[m][n],double y[m][n],
                     double x_ih[m][n-1],double y_ih[m][n-1],
                     double x_jh[m-1][n],double y_jh[m-1][n]);
void get_u_v_potential_fullp(int m,int n,double phi[m][n],double x[m][n],
                             double y[m][n],double J[m][n],double A2[m][n],
                             double A3[m][n],double u[m][n],double v[m][n],
                             double q[m][n],int op);
void initialize_fullp(int m,int n,double phi[m][n],double x[m][n],
                      double y[m][n],fullp_prmtrs *fp_prmtrs);
void L_phi_fullp(int m,int n,double L_phi[m][n],double L_phi_terms_ih[m][n-1],
                 double L_phi_terms_jh[m][n-1]);
void L_phi_fullp_der_terms_ih(int m,int n,double L_phi_terms_ih[m][n-1],
                              double rho_til[m][n],double contraU_ih[m][n],
                              double J_ih[m][n]);
void L_phi_fullp_der_terms_jh(int m,int n,double L_phi_terms_jh[m][n-1],
                              double rho_bar[m][n],double contraV_jh[m][n],
                              double J_jh[m][n]);
// void local_q_dervs(double *dphi_dksi,double *dphi_deta,int m,int n,
//                    double phi[m][n],double A2[m][n],double A3[m][n],
//                    int i,int j,int af2_flip);
double max2(double n1,double n2);
// double mean2(double n1,double n2);
double mean4_j(int m,int n,double rho[m][n-1],int i,int j);
// void metric_terms(double *ksix,double *ksiy,double *etax,double *etay,int m,
//                   int n,double J[m][n],double x[m][n],double y[m][n],
//                   int i,int j);
double norm_L2(int m,int n,double L_phi[m][n]);
void nu_switch(int m,int n,double nu[m][n],double contraUV[m][n],
               double rho[m][n],double C,int axis);
double rho_coeffs(int m,int n,double rho_tb[m][n],double rho[m][n],
                  double contraUV[m][n],double nu[m][n],int op,int axis);
int rs_idx(double contraUV);
void three_point_pol2_extrp(int m,int n,double A[m][n],double x[m][n],
                            double y[m][n],int i,int end);
void two_point_linear_extrp(int m,int n,double A[m][n],double x[m][n],
                            double y[m][n],int i,int end);
// void solve_adi_2d_rectangular_fullp(int m,int n,double phi[m][n],double J[m][n],
//                                     double A1[m][n],double A2[m][n],
//                                     double A3[m][n],sim_prmtrs *config,
//                                     fullp_prmtrs *fp_prmtrs);
void solve_af2_2d_rectangular_fullp(int m,int n,double phi[m][n],double J[m][n],
                                    double A1[m][n],double A2[m][n],
                                    double A3[m][n],double J_ih[m][n],
                                    double A1_ih[m][n],double A2_ih[m][n],
                                    double A3_ih[m][n],double J_jh[m][n],
                                    double A1_jh[m][n],double A2_jh[m][n],
                                    double A3_jh[m][n],sim_prmtrs *config,
                                    fullp_prmtrs *fp_prmtrs);

#endif