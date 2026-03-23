#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/mesh.h"
#include"../include/num_methods.h"
#include"../include/strings.h"
#include"./bi_air_lib.h"

/*
- gigiaero, 22/03/2026, 2131 hours
*/
void biconvex_airfoil_mesh(bi_air_mesh *b_a_m,double *x, double*y){

  double delta_x = 1./(b_a_m->ITE -b_a_m->ILE);

  for(int i=b_a_m->ILE-1;i<b_a_m->ITE;i++)
    x[i] = (i - b_a_m->ILE)*delta_x;

  for(int i=b_a_m->ITE;i<b_a_m->IMAX;i++)
    x[i] = x[i-1] + (x[i-1] - x[i-2])*b_a_m->XSF;

  for(int i=b_a_m->ILE-2;i>=0;i--)
    x[i] = x[i+1] + (x[i+1] - x[i+2])*b_a_m->XSF;

  y[0] = -delta_x/2.;
  y[1] = -y[0];

  for(int j=2;j<b_a_m->JMAX;j++)
    y[j] = y[j-1] + (y[j-1] - y[j-2])*b_a_m->YSF;
}

int main(){

  // // Solution configurations
  // sim_parameters config;
  // config.Ntype = 1;
  // config.max_iter = 5000;
  // config.qtimes = 10;
  // config.save_last_only = 1;
  // config.eps = 1.e-5; // Convergence criterion
  // char output_file[] = "results/laplace2d";

  // Mesh
  bi_air_mesh b_a_mesh;
  b_a_mesh.ILE = 11;   // Leading edge
  b_a_mesh.ITE = 31;   // Trailing edge
  b_a_mesh.IMAX = 41;  // Number of points along x
  b_a_mesh.JMAX = 12;  // Number of points along y
  b_a_mesh.XSF = 1.25; // Stretching factor, x
  b_a_mesh.YSF = 1.25; // Stretching factor, y

  // Physical properties
  // double alpha = 1.; // CONFIRMAR ISTO
  // double lx = 4.;
  // double ly = 2.;

  // // Boundary conditions
  // int num_b_cs = 4;
  // b_conditions_2d b_c[num_b_cs];
  // // Down
  // b_c[0].type = 'D'; 
  // b_c[0].val = 0.;
  // b_c[0].axis = 1;
  // b_c[0].position = 0;
  // b_c[0].range[0] = 1;
  // b_c[0].range[1] = n-1;
  // // Right
  // b_c[1].type = 'D'; 
  // b_c[1].val = -10.;
  // b_c[1].axis = 2;
  // b_c[1].position = n-1;
  // b_c[1].range[0] = 1;
  // b_c[1].range[1] = m-1;
  // // Up
  // b_c[2].type = 'D'; 
  // b_c[2].val = 50.;
  // b_c[2].axis = 1;
  // b_c[2].position = m-1;
  // b_c[2].range[0] = 1;
  // b_c[2].range[1] = n-1;
  // // Left
  // b_c[3].type = 'D'; 
  // b_c[3].val = 10.;
  // b_c[3].axis = 2;
  // b_c[3].position = 0;
  // b_c[3].range[0] = 1;
  // b_c[3].range[1] = m-1;

  // // Initialize mesh and temperature array
  double *x = malloc(sizeof(double)*b_a_mesh.IMAX);
  double *y = malloc(sizeof(double)*b_a_mesh.JMAX);
  // double (*T)[n] = calloc(m,sizeof *T);
  // double delta_x = lx/(n - 1.);
  // double delta_y = ly/(m - 1.);
  // zeros_2d_array(m,n,T);
  biconvex_airfoil_mesh(&b_a_mesh,x,y);

  print_1d_array(b_a_mesh.IMAX,x);
  putchar('\n');
  print_1d_array(b_a_mesh.JMAX,y);

  

  // // Initialize boundary conditions
  // apply_b_cs(m,n,T,num_b_cs,b_c);

  // // Solve
  // config.casename = malloc(sizeof(char)*200);
  // strcat(config.casename,output_file);
  // evaluate_delta_form(m,n,T,x,y,&config);

  free(x);
  free(y);
  // free(T);
  // free(config.casename);

  return 0;
}