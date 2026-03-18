#include<stdio.h>
#include<stdlib.h>
#include"../../include/1d_arrays.h"
#include"../../include/2d_arrays.h"

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

  // print_1d_array(n-1,a);
  // print_1d_array(n,b);
  // print_1d_array(n-1,c);
  // print_1d_array(n-2,d);

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

  putchar('k');
  print_1d_array(n-1,k);
  putchar('w');
  print_1d_array(n,w);
  putchar('x');
  print_1d_array(n-1,x);
  putchar('z');
  print_1d_array(n-2,z);

  // Solve for y
  y[0] = f[0]/w[0];
  for(int i=1;i<n;i++)
    y[i] = (f[i] - k[i-1]*y[i-1])/w[i];

  putchar('\n');
  print_1d_array(n,y);

  // Solve for u
  u[n-1] = y[n-1];
  u[n-2] = y[n-2] - x[n-2]*u[n-1];
  // printf("%f | %f | %f\n",y[n-2],x[n-2],u[n-1]);
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

int main(){

  int n = 10;
  double A[10][10] = {{3.,3.,6.,0.,0.,0.,0.,0.,0.,0.}
  ,{8.,9.,6.,2.,0.,0.,0.,0.,0.,0.}
  ,{0.,1.,4.,2.,4.,0.,0.,0.,0.,0.}
  ,{0.,0.,2.,7.,9.,3.,0.,0.,0.,0.}
  ,{0.,0.,0.,5.,1.,6.,3.,0.,0.,0.}
  ,{0.,0.,0.,0.,7.,8.,4.,3.,0.,0.}
  ,{0.,0.,0.,0.,0.,3.,1.,1.,1.,0.}
  ,{0.,0.,0.,0.,0.,0.,9.,1.,3.,2.}
  ,{0.,0.,0.,0.,0.,0.,0.,2.,7.,4.}
  ,{0.,0.,0.,0.,0.,0.,0.,0.,9.,9.}};
  double f[] = {2.,3.,4.,8.,1.,3.,8.,6.,7.,3.};
  
  
  
  
  double u[n];

  diagonal_matrix_solver(n,A,f,u);

  print_1d_array(n,u);

  return 0;

}


// [[9. 1. 1. 0. 0. 0. 0. 0. 0. 0.]
//  [1. 7. 6. 9. 0. 0. 0. 0. 0. 0.]
//  [0. 5. 3. 5. 6. 0. 0. 0. 0. 0.]
//  [0. 0. 3. 7. 7. 6. 0. 0. 0. 0.]
//  [0. 0. 0. 2. 2. 8. 2. 0. 0. 0.]
//  [0. 0. 0. 0. 1. 1. 8. 8. 0. 0.]
//  [0. 0. 0. 0. 0. 9. 3. 5. 7. 0.]
//  [0. 0. 0. 0. 0. 0. 6. 5. 2. 1.]
//  [0. 0. 0. 0. 0. 0. 0. 6. 3. 2.]
//  [0. 0. 0. 0. 0. 0. 0. 0. 4. 4.]]
// [7. 2. 7. 7. 7. 2. 4. 6. 7. 2.]
// [[9. 1. 1. 0. 0. 0. 0. 0. 0. 0.]
//  [1. 7. 6. 9. 0. 0. 0. 0. 0. 0.]
//  [0. 5. 3. 5. 6. 0. 0. 0. 0. 0.]
//  [0. 0. 3. 7. 7. 6. 0. 0. 0. 0.]
//  [0. 0. 0. 2. 2. 8. 2. 0. 0. 0.]
//  [0. 0. 0. 0. 1. 1. 8. 8. 0. 0.]
//  [0. 0. 0. 0. 0. 9. 3. 5. 7. 0.]
//  [0. 0. 0. 0. 0. 0. 6. 5. 2. 1.]
//  [0. 0. 0. 0. 0. 0. 0. 6. 3. 2.]
//  [0. 0. 0. 0. 0. 0. 0. 0. 4. 4.]]
// True

// [ 2.81662344e+00  2.26953748e+00 -2.06191484e+01  1.18901672e+01
//  -3.23513016e-01 -2.01818900e+00  6.10182178e-03  5.36610931e-01
//   2.78033442e+00 -2.28033442e+00]
// [ 2.81662344e+00  2.26953748e+00 -2.06191484e+01  1.18901672e+01
//  -3.23513016e-01 -2.01818900e+00  6.10182178e-03  5.36610931e-01
//   2.78033442e+00 -2.28033442e+00]
