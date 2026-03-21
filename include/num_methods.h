#ifndef _lib_num_methods_
#define _lib_num_methods_

/*
Defined in the context of laplace2d
-gigiaero, 20/03/2026, 0059 hours
*/
typedef struct s_parameters{
  int Ntype;
  int max_iter;
  int qtimes;
  double eps;
  char casename[200];
}sim_parameters;

double delta_xy(double *xy,int i);
void diagonal_matrix_solver(int n,double A[n][n],double *f,double *u);
double scheme_der2_o2_central(double phi_ip1,double phi_i,double phi_im1,
                              double x_ip1,double x_i,double x_im1);

#endif