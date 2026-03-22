#ifndef _lib_num_methods_
#define _lib_num_methods_

/*
Defined in the context of laplace2d
- gigiaero, 20/03/2026, 0059 hours
*/
typedef struct s_parameters{
  int Ntype;
  int max_iter;
  int qtimes; // Set to 0 to save only the final result
  int save_last_only;
  double eps;
  char *casename;
}sim_parameters;

void calc_residual(double *phi_elem,double *phi_old_elem,double *N_elem,
                  double *res_elem);
int check_num_digits_int(int *num);
double delta_xy(double *xy,int i);
void diagonal_matrix_solver(int n,double A[n][n],double *f,double *u);
void evaluate_delta_form(int m,int n,double phi[m][n],double *x,double *y,
                         sim_parameters *config);
void N_p_jacobi(double *N,double *x,double *y,int i,int j);
// void save_results_qtimes(int m,int n,double phi[m][n],int *iter,int *qtimes,
//                          char *buffer,char *filename_save,int *str_end_idx);
void save_results_qtimes(int m,int n,double phi[m][n],int *iter,sim_parameters *s_p,
  char *buffer,char *filename_save,int *str_end_idx);
double scheme_der2_o2_central(double phi_ip1,double phi_i,double phi_im1,
                              double x_ip1,double x_i,double x_im1);
void scheme_der2_o2_central_var_deltas_xy(double *f,int m,int n,
                                          double phi[m][n],double *x,double *y,
                                          int i,int j);
void solve_p_jacobi_2d_rectangular(int m,int n,double phi[m][n],double *x,
                                   double *y,sim_parameters *config);

#endif