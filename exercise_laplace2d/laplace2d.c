#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/mesh.h"
#include"../include/num_methods.h"
#include"../include/strings.h"


int main(){

  // Solution configurations
  sim_prmtrs config;
  config.Ntype = 4;
  config.r = 1;
  config.max_iter = 50000;
  config.qtimes = 10;
  config.save_i_c = 0;
  config.save_last_only = 1;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/laplace2d";

  // Mesh
  int m = 50;
  int n = 50;

  // Physical properties
  // double alpha = 1.;
  double lx = 4.;
  double ly = 1.;

  // Boundary conditions
  int num_b_c = 4;
  b_conditions_2d b_c[num_b_c];
  for(int i=0;i<num_b_c;i++)
    b_c[i].val = malloc(sizeof(double));
  int counter = 0;
  // Down
  b_c[counter].type = 'D'; 
  b_c[counter].val[0] = 0.;
  b_c[counter].axis = 1;
  b_c[counter].position = 0;
  b_c[counter].range[0] = 1;
  b_c[counter].range[1] = n-2;
  counter++;
  // Right
  b_c[counter].type = 'D'; 
  b_c[counter].val[0] = -100.;
  b_c[counter].axis = 2;
  b_c[counter].position = n-1;
  b_c[counter].range[0] = 1;
  b_c[counter].range[1] = m-2;
  counter++;
  // Up
  b_c[counter].type = 'D'; 
  b_c[counter].val[0] = 0.;
  b_c[counter].axis = 1;
  b_c[counter].position = m-1;
  b_c[counter].range[0] = 1;
  b_c[counter].range[1] = n-2;
  counter++;
  // Left
  b_c[counter].type = 'D'; 
  b_c[counter].val[0] = 100.;
  b_c[counter].axis = 2;
  b_c[counter].position = 0;
  b_c[counter].range[0] = 1;
  b_c[counter].range[1] = m-2;

  // Reapplied boundary conditions
  int num_b_c_re = 0;
  b_conditions_2d b_c_re[num_b_c_re];
  for(int i=0;i<num_b_c_re;i++)
    b_c_re[i].val = malloc(sizeof(double));
  counter = 0;
  // // Down
  // b_c_re[counter].type = 'N'; 
  // b_c_re[counter].val[0] = -10.;
  // b_c_re[counter].axis = 1;
  // b_c_re[counter].position = 0;
  // b_c_re[counter].range[0] = 1;
  // b_c_re[counter].range[1] = n-2;
  // counter++;
  // // Right
  // b_c_re[counter].type = 
  // b_c_re[counter].val[0] = 0.;
  // b_c_re[counter].axis = 2;
  // b_c_re[counter].position = n-1;
  // b_c_re[counter].range[0] = 1;
  // b_c_re[counter].range[1] = m-2;
  // counter++;
  // Up
  // b_c_re[counter].type = 'N';
  // b_c_re[counter].val[0] = -10.;
  // b_c_re[counter].axis = 1;
  // b_c_re[counter].position = m-1;
  // b_c_re[counter].range[0] = 1;
  // b_c_re[counter].range[1] = n-2;
  // counter++;
  // // Left
  // b_c_re[counter].type = 
  // b_c_re[counter].val[0] = 0.;
  // b_c_re[counter].axis = 2;
  // b_c_re[counter].position = 0;
  // b_c_re[counter].range[0] = 1;
  // b_c_re[counter].range[1] = m-2;

  // Initialize mesh and temperature array
  double *x = malloc(sizeof(double)*n);
  double *y = malloc(sizeof(double)*m);
  double (*T)[n] = calloc(m,sizeof *T);
  double delta_x = lx/(n - 1.);
  double delta_y = ly/(m - 1.);
  zeros_2d_array(m,n,T);
  uniform_rectangular_mesh(m,n,delta_x,delta_y,x,y);

  // Initialize boundary conditions
  apply_b_c(m,n,T,num_b_c,b_c,x,y);

  // Solve
  config.casename = malloc(sizeof(char)*200);
  strcpy(config.casename,output_file);
  evaluate_delta_form(m,n,T,x,y,&config,num_b_c_re,b_c_re);

  for(int i=0;i<num_b_c;i++)
    free(b_c[i].val);
  for(int i=0;i<num_b_c_re;i++)
    free(b_c_re[i].val);
  free(x);
  free(y);
  free(T);
  free(config.casename);

  return 0;
}