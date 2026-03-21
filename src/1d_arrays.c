#include<stdio.h>
#include"../include/1d_arrays.h"

/*
- gigiaero, 12/03/2026, 2133 hours
*/
void fill_1d_array(int n,double *v){
  int k = 0;
  for(int i=0;i<n;i++){
    v[i] = (double) k;
    k++;
  }
}

/*
- gigiaero, 12/03/2026, 1704 hours
*/
void ones_1d_array(int n,double *v){
  for(int i=0;i<n;i++)
    v[i] = 1.;
}

/*
- gigiaero, 12/03/2026, 1708 hours
*/
void print_1d_array(int n,double *v){
  for(int i=0;i<n;i++)
    printf("%f ",v[i]);

  putchar('\n');
}

void print_1d_array_to_file(int n,double *v,char *filename){
  FILE *output;

  output = fopen(filename,"w");

  for(int i=0;i<n;i++)
    fprintf(output,"%f ",v[i]);

  fclose(output);
}

/*
- gigiaero, 12/03/2026, 1704 hours
*/
void zeros_1d_array(int n,double *v){
  for(int i=0;i<n;i++)
    v[i] = 0.;
}