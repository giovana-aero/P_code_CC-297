#include<stdio.h>
#include"../include/2d_arrays.h"

/*
Copy array A to B
- gigiaero, 19/03/2026, 2200 hours
*/
void copy_2d_array(int m,int n,double A[m][n],double B[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      B[j][i] = A[j][i];
    
  }
}

/* 
Fills a 2d array with increasing values
- gigiaero, 10/03/2026, 1344 hours
*/
void fill_2d_array(int m,int n,double A[m][n]){
  int k = 0;
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++){
      A[j][i] = (double) k;
      k++;
    }
  }
}

/*
- gigiaero, 10/03/2026, 1345 hours
*/
void ones_2d_array(int m,int n,double A[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      A[j][i] = 1.;
  }
}

/*
- gigiaero, 10/03/2026, 1346 hours
*/
void print_2d_array(int m,int n,double A[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      printf("%f ",A[j][i]);
    putchar('\n');
  }
}

// print to file goes here

/*
- gigiaero, 10/03/2026, 1346 hours
*/
void zeros_2d_array(int m,int n,double A[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      A[j][i] = 0.;
  }
}
