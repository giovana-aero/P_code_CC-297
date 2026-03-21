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
  // if(invert_y){
  //   for(int j=m-1;j>=0;j--){
  //     for(int i=0;i<n;i++)
  //       printf("%f ",A[j][i]);
  //     putchar('\n');
  //   }
  // }
  // else{
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      printf("%f ",A[j][i]);
    putchar('\n');
  }
  // }
}

/*
- gigiaero, 20/03/2026, 0058 hours
*/
void print_2d_array_to_file(int m,int n,double A[m][n],char *filename,
                            int invert_y){
  FILE *output;

  output = fopen(filename,"w");

  if(invert_y){
    for(int j=m-1;j>=0;j--){
      for(int i=0;i<n;i++)
        fprintf(output,"%.6f ",A[j][i]);
      
      fprintf(output,"\n");
    }  
  }
  else{
    for(int j=0;j<m;j++){
      for(int i=0;i<n;i++)
        fprintf(output,"%.6f ",A[j][i]);
      
      fprintf(output,"\n");
    }
  }

  fclose(output);
}

/*
- gigiaero, 10/03/2026, 1346 hours
*/
void zeros_2d_array(int m,int n,double A[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      A[j][i] = 0.;
  }
}
