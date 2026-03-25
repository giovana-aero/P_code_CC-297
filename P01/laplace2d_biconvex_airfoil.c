#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/bi_air_lib.h"
#include"../include/mesh.h"
#include"../include/num_methods.h"
#include"../include/strings.h"

int main(){

  // Solution configurations
  sim_parameters config;
  config.Ntype = 1;
  config.max_iter = 50000;
  config.qtimes = 10;
  config.save_last_only = 1;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/bi_air";

  // Mesh
  bi_air_phys_mesh b_a_mesh;
  b_a_mesh.ILE = 11;   // Leading edge
  b_a_mesh.ITE = 31;   // Trailing edge
  b_a_mesh.IMAX = 41;  // Number of points along x
  b_a_mesh.JMAX = 12;  // Number of points along y
  b_a_mesh.XSF = 1.25; // Stretching factor, x
  b_a_mesh.YSF = 1.25; // Stretching factor, y

  // Physical properties
  b_a_mesh.uinf = 1.;
  b_a_mesh.t = 0.05;

  // Initialize mesh and velocity potential array
  double *x = malloc(sizeof(double)*b_a_mesh.IMAX);
  double *y = malloc(sizeof(double)*b_a_mesh.JMAX);
  double (*phi)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *phi);
  double (*u)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *u);
  double (*v)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *v);
  biconvex_airfoil_mesh(&b_a_mesh,x,y);

  // print_1d_array(b_a_mesh.IMAX,x);
  // putchar('\n');
  // print_1d_array(b_a_mesh.JMAX,y);


  // Boundary conditions
  int num_b_c = 3;
  b_conditions_2d b_c[num_b_c];
  int counter = 0;
  // // Down, airfoil                                   [CORRECTION PENDING]
  // b_c[counter].type = 'd'; 
  // b_c[counter].axis = 1;
  // b_c[counter].position = 0;
  // b_c[counter].range[0] = b_a_mesh.ILE - 1;
  // b_c[counter].range[1] = b_a_mesh.ITE - 1;
  // b_c[counter].val = malloc(sizeof(double)*(b_c[counter].range[1] - b_c[counter].range[0] + 1));
  // bi_air_dirichlet_vals_wall(x,&b_c[0],uinf,t);
  // counter++;
  // // Down, upstream of the leading edge                   [CORRECTION PENDING]
  // b_c[counter].type = 'd'; 
  // b_c[counter].axis = 1;
  // b_c[counter].position = 0;
  // b_c[counter].range[0] = 0;
  // b_c[counter].range[1] = b_a_mesh.ILE - 2;
  // b_c[counter].val = malloc(sizeof(double)*(b_c[counter].range[1] - b_c[counter].range[0] + 1));
  // bi_air_dirichlet_vals_free(x,&b_c[1],uinf);
  // counter++;
  // // Down, downstream of the leading edge                  [CORRECTION PENDING]
  // b_c[counter].type = 'd'; 
  // b_c[counter].axis = 1;
  // b_c[counter].position = 0;
  // b_c[counter].range[0] = b_a_mesh.ITE;
  // b_c[counter].range[1] = b_a_mesh.IMAX - 1;
  // b_c[counter].val = malloc(sizeof(double)*(b_c[counter].range[1] - b_c[counter].range[0] + 1));
  // bi_air_dirichlet_vals_free(x,&b_c[2],uinf);
  // counter++;
  // Left (inlet)
  b_c[counter].type = 'D';
  b_c[counter].axis = 2;
  b_c[counter].position = 0;
  b_c[counter].range[0] = 0;
  b_c[counter].range[1] = b_a_mesh.JMAX - 1;
  b_c[counter].val = malloc(sizeof(double));
  b_c[counter].val[0] = b_a_mesh.uinf*x[0];
  counter++;
  // Right (outflow)
  b_c[counter].type = 'D'; 
  b_c[counter].axis = 2;
  b_c[counter].position = b_a_mesh.IMAX - 1;
  b_c[counter].range[0] = 0;
  b_c[counter].range[1] = b_a_mesh.JMAX - 1;
  b_c[counter].val = malloc(sizeof(double));
  b_c[counter].val[0] = b_a_mesh.uinf*x[b_a_mesh.IMAX-1];
  counter++;
  // Up (open boundary)
  b_c[counter].type = 'd'; 
  b_c[counter].axis = 1;
  b_c[counter].position = b_a_mesh.JMAX - 1;
  b_c[counter].range[0] = 1;
  b_c[counter].range[1] = b_a_mesh.IMAX - 2;
  b_c[counter].val = malloc(sizeof(double)*(b_c[counter].range[1] - b_c[counter].range[0] + 1));
  bi_air_dirichlet_vals_free(x,&b_c[counter],b_a_mesh.uinf);
  
  // Reapplied boundary conditions
  int num_b_c_r = 0;
  b_conditions_2d b_c_r[num_b_c_r];
  // // Down, airfoil
  // b_c[0].type = 'd'; 
  // b_c[0].axis = 1;
  // b_c[0].position = 0;
  // b_c[0].range[0] = b_a_mesh.ILE - 1;
  // b_c[0].range[1] = b_a_mesh.ITE - 1;
  // b_c[0].val = malloc(sizeof(double)*(b_c[0].range[1] - b_c[0].range[0] + 1));
  // bi_air_dirichlet_vals_wall(x,&b_c[0],uinf,t);
  // // Down, upstream of the leading edge
  // b_c[1].type = 'd'; 
  // b_c[1].axis = 1;
  // b_c[1].position = 0;
  // b_c[1].range[0] = 0;
  // b_c[1].range[1] = b_a_mesh.ILE - 2;
  // b_c[1].val = malloc(sizeof(double)*(b_c[1].range[1] - b_c[1].range[0] + 1));
  // bi_air_dirichlet_vals_free(x,&b_c[1],uinf);
  // // Down, downstream of the leading edge
  // b_c[2].type = 'd'; 
  // b_c[2].axis = 1;
  // b_c[2].position = 0;
  // b_c[2].range[0] = b_a_mesh.ITE;
  // b_c[2].range[1] = b_a_mesh.IMAX - 1;
  // b_c[2].val = malloc(sizeof(double)*(b_c[2].range[1] - b_c[2].range[0] + 1));
  // bi_air_dirichlet_vals_free(x,&b_c[2],uinf);

  // Initialize boundary conditions
  apply_b_c(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,num_b_c,b_c,x,y);

  // Initialize initial conditions
  for(int j=1;j<b_a_mesh.JMAX-1;j++)
    copy_1d_array_range(1,b_a_mesh.IMAX-1,phi[b_a_mesh.JMAX-1],phi[j]);

  // char filename[] = "test.txt";
  // print_2d_array_to_file(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,filename,1);

  // Solve
  config.casename = malloc(sizeof(char)*200);
  strcpy(config.casename,output_file);
  evaluate_delta_form_bi_air(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,x,y,&config,
                             &b_a_mesh);

  // save_mesh(b_a_mesh.JMAX,b_a_mesh.IMAX,x,y,config.casename);




  for(int i=0;i<num_b_c;i++)
    free(b_c[i].val);
  free(x);
  free(y);
  free(phi);
  free(config.casename);

  return 0;
}