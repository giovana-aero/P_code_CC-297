#include<stdio.h>
#include<stdlib.h>

int main(){

  int n = 10;

  double a[] = {1.,9.,5.,2.,0.,2.,6.,4.,0.,8.};
  double b[] = {5.,3.,8.,6.,7.,1.,1.,1.,1.,5.};
  double c[] = {5.,2.,6.,5.,3.,9.,5.,8.,4.,4.};
  double f[] = {8.,9.,3.,3.,5.,6.,5.,0.,3.,2.};
  double u[n];

  tridiagonal_pmatrix_solver(n,a,b,c,f,u);

  for(int i=0;i<n;i++)
    printf("%f ",u[i]);
  putchar('\n');

  return 0;
}