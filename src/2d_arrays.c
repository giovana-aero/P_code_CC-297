#include<stdio.h>
#include"../include/2d_arrays.h"

/* 
Fills a 2d array with increasing values
- gigiaero, 10/03/2026, 1344 hours
*/
void fill_2d_array(int m,int n,double A[m][n]){
  int k = 0;
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      A[i][j] = (double) k;
      k++;
    }
  }
}

/*
- gigiaero, 10/03/2026, 1345 hours
*/
void ones_2d(int m,int n,double A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      A[i][j] = 1.;
  }
}

/*
- gigiaero, 10/03/2026, 1346 hours
*/
void print_2d_array(int m,int n,double A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      printf("%f ",A[i][j]);
    putchar('\n');
  }
}

/*
- gigiaero, 10/03/2026, 1346 hours
*/
void zeros_2d(int m,int n,double A[m][n]){
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++)
      A[i][j] = 0.;
  }
}
