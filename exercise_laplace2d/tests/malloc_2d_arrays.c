#include<stdio.h>
#include<stdlib.h>
#include"../../include/2d_arrays.h"

void execute(int m,int n){
  double **A;
  A = (double**)malloc(sizeof(double*)*n + sizeof(double)*m*n);

  // int k = 0;
  // for(int j=0;j<m;j++){
  //   for(int i=0;i<n;i++){
  //     A[j][i] = (double) k;
  //     k++;
  //   }
  // }

  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      printf("%f ",A[j][i]);
    putchar('\n');
  }

  free(A);
}

int main(){

  int m = 4;
  int n = 6;

  execute(m,n);

  return 0;
}