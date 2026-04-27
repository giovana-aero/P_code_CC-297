#include<stdio.h>
#include<stdlib.h>
#include"../../include/num_methods.h"

int main(){
  int n = 15;

  double a[] = {5.,8.,1.,6.,3.,8.,9.,7.,5.,1.,4.,8.,3.,3.};
  double b[] = {3.,2.,9.,4.,9.,8.,5.,1.,4.,7.,2.,7.,2.,9.,6.};
  double c[] = {8.,6.,8.,3.,6.,5.,4.,3.,4.,1.,2.,2.,3.,1.};
  double f[] = {8.,9.,1.,3.,3.,7.,8.,3.,5.,8.,1.,5.,8.,2.,2.};
  double u[n];

  tridiagonal_matrix_solver(n,a,b,c,f,u);

  for(int i=0;i<n;i++)
    printf("%f ",u[i]);
  putchar('\n');

  return 0;
}