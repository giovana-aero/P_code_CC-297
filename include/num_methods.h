#ifndef _lib_num_methods_
#define _lib_num_methods_

void diagonal_matrix_solver(int n,double A[n][n],double *f,double *u);
double scheme_der2_o2_central(double phi_ip1,double phi_i,double phi_im1,
                              double x_ip1,double x_i,double x_im1);

#endif