#include<stdlib.h>
#include<string.h>
#include"../include/1d_arrays.h"
#include"../include/mesh.h"
#include"../include/strings.h"

/*
- gigiaero, 22/03/2026, 1034 hours
*/
void save_mesh(int m,int n,double *x,double *y,char *casename){
  char *filename = malloc(sizeof(char)*200);
  int str_end_idx;

  find_str_end(filename,&str_end_idx);

  strcpy(filename,casename);
  find_str_end(filename,&str_end_idx);
  strcat(filename,"_mesh_x.msh");
  print_1d_array_to_file(n,x,filename);
  filename[str_end_idx] = '\0';
  strcat(filename,"_mesh_y.msh");
  print_1d_array_to_file(m,y,filename);

  free(filename);
}

void uniform_rectangular_mesh(int m,int n,double delta_x,double delta_y,
  double *x,double *y){
  x[0] = 0.;
  y[0] = 0.;

  for(int i=1;i<n;i++)
  x[i] = x[i-1] + delta_x;

  for(int j=1;j<m;j++)
  y[j] = y[j-1] + delta_y;
}