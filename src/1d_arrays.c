#include<stdio.h>
#include"../include/1d_arrays.h"

/*
Copies array v1 to v2 (of same size), with ranges specified with an initial 
index and a final index
- gigiaero, 23/03/2026, 1218 hours
*/
void copy_1d_array_range(int i_1,int i_f,double *v1,double *v2){
  for(int i=i_1;i<=i_f;i++)
    v2[i] = v1[i];
}

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
    printf("%e ",v[i]);

  putchar('\n');
}

/*
- gigiaero, 22/03/2026, 0950 hours
*/
void print_1d_array_int(int n,int *v){
  for(int i=0;i<n;i++)
    printf("%d ",v[i]);

  putchar('\n');
}

/*
- gigiaero, 20/03/2026, 2300 hours
*/
void print_1d_array_to_file(int n,double *v,char *filename){
  FILE *output;

  output = fopen(filename,"w");

  for(int i=0;i<n;i++)
    fprintf(output,"%.6f\n",v[i]);

  fclose(output);
}

/*
- gigiaero, 12/03/2026, 1704 hours
*/
void zeros_1d_array(int n,double *v){
  for(int i=0;i<n;i++)
    v[i] = 0.;
}