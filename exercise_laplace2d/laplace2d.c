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
  sim_parameters config;
  config.Ntype = 1;
  config.max_iter = 10000;
  config.qtimes = 10;
  config.save_last_only = 1;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/laplace2d";

  // Mesh
  int m = 100;
  int n = 100;

  // Physical properties
  // double alpha = 1.; // CONFIRMAR ISTO
  double lx = 4.;
  double ly = 2.;

  // Boundary conditions
  int num_b_cs = 4;
  b_conditions_2d b_c[num_b_cs];
  // Down
  b_c[0].type = 'D'; 
  b_c[0].val = 0.;
  b_c[0].axis = 1;
  b_c[0].position = 0;
  b_c[0].range[0] = 1;
  b_c[0].range[1] = n-1;
  // Right
  b_c[1].type = 'D'; 
  b_c[1].val = -10.;
  b_c[1].axis = 2;
  b_c[1].position = n-1;
  b_c[1].range[0] = 1;
  b_c[1].range[1] = m-1;
  // Up
  b_c[2].type = 'D'; 
  b_c[2].val = 50.;
  b_c[2].axis = 1;
  b_c[2].position = m-1;
  b_c[2].range[0] = 1;
  b_c[2].range[1] = n-1;
  // Left
  b_c[3].type = 'D'; 
  b_c[3].val = 10.;
  b_c[3].axis = 2;
  b_c[3].position = 0;
  b_c[3].range[0] = 1;
  b_c[3].range[1] = m-1;

  // Initialize mesh and temperature array
  double *x = malloc(sizeof(double)*n);
  double *y = malloc(sizeof(double)*m);
  double (*T)[n] = calloc(m,sizeof *T);
  double delta_x = lx/(n - 1.);
  double delta_y = ly/(m - 1.);
  zeros_2d_array(m,n,T);
  uniform_rectangular_mesh(m,n,delta_x,delta_y,x,y);

  // Initialize boundary conditions
  apply_b_cs(m,n,T,num_b_cs,b_c);

  // Solve
  config.casename = malloc(sizeof(char)*200);
  strcat(config.casename,output_file);
  evaluate_delta_form(m,n,T,x,y,&config);

  free(x);
  free(y);
  free(T);
  free(config.casename);

  return 0;
}