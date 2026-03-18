#include<stdio.h>
#include<stdlib.h>
#include"../include/num_methods.h"

/*
This was based on series 02's question 5. Also works if some of the diagonals 
are zero (except the main one)
- gigiaero, 18/03/2026, 1620 hours
*/
void diagonal_matrix_solver(int n,double A[n][n],double *f,double *u){

  // Original matrix
  double *a = (double *)malloc(sizeof(double)*(n-1));
  double *b = (double *)malloc(sizeof(double)*n);
  double *c = (double *)malloc(sizeof(double)*(n-1));
  double *d = (double *)malloc(sizeof(double)*(n-2));
  // Decomposed matrices
  double *k = (double *)malloc(sizeof(double)*(n-1));
  double *w = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*(n-1));
  double *z = (double *)malloc(sizeof(double)*(n-2));
  // Intermediary solution
  double *y = (double *)malloc(sizeof(double)*n);
  

  if(a == NULL){
    puts("diagonal_matrix_solver: Allocation failed (a)");
    exit(1);
  }
  if(b == NULL){
    puts("diagonal_matrix_solver: Allocation failed (b)");
    exit(1);
  }
  if(c == NULL){
    puts("diagonal_matrix_solver: Allocation failed (c)");
    exit(1);
  }
  if(d == NULL){
    puts("diagonal_matrix_solver: Allocation failed (d)");
    exit(1);
  }
  if(k == NULL){
    puts("diagonal_matrix_solver: Allocation failed (k)");
    exit(1);
  }
  if(w == NULL){
    puts("diagonal_matrix_solver: Allocation failed (w)");
    exit(1);
  }
  if(x == NULL){
    puts("diagonal_matrix_solver: Allocation failed (x)");
    exit(1);
  }
  if(z == NULL){
    puts("diagonal_matrix_solver: Allocation failed (z)");
    exit(1);
  }
  if(y == NULL){
    puts("diagonal_matrix_solver: Allocation failed (y)");
    exit(1);
  }

  // Get coefficients of the original matrix
  for(int i=0;i<n;i++){
    b[i] = A[i][i];

    if(i < n-1)
      a[i] = A[i+1][i];

    if(i > 0)
      c[i-1] = A[i-1][i];

    if(i > 1)
      d[i-2] = A[i-2][i];
  }

  // Calculate coefficients of decomposed matrices
  k[0] = a[0];
  w[0] = b[0];
  x[0] = c[0]/w[0];
  z[0] = d[0]/w[0];
  for(int i=1;i<n;i++){
    if(i < (n - 1))
      k[i] = a[i];

    w[i] = b[i] - k[i-1]*x[i-1];

    if(i < (n - 1))
      x[i] = (c[i] - k[i-1]*z[i-1])/w[i];
    
    if(i < (n - 2))
      z[i] = d[i]/w[i];
  }

  // Solve for y
  y[0] = f[0]/w[0];
  for(int i=1;i<n;i++)
    y[i] = (f[i] - k[i-1]*y[i-1])/w[i];

  // Solve for u
  u[n-1] = y[n-1];
  u[n-2] = y[n-2] - x[n-2]*u[n-1];
  for(int i=n-3;i>=0;i--)
    u[i] = y[i] - x[i]*u[i+1] - z[i]*u[i+2];

  free(a);
  free(b);
  free(c);
  free(d);
  free(k);
  free(w);
  free(x);
  free(z);
  free(y);
}

/*
Scheme: second derivative, second order, central
- gigiaero, 13/03/2026, 2315 hours
*/
double scheme_der2_o2_central(double phi_ip1,double phi_i,double phi_im1,
                              double x_ip1,double x_i,double x_im1){

  return 2./(x_ip1 - x_im1)*((phi_ip1 - phi_i)/(x_ip1 - x_i) - 
         (phi_i - phi_im1)/(x_i - x_im1));

}