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

// /*
// - gigiaero, 22/03/2026, 2131 hours
// */
// void biconvex_airfoil_mesh(bi_air_mesh *b_a_m,double *x, double*y){

//   double delta_x = 1./(b_a_m->ITE -b_a_m->ILE);

//   for(int i=b_a_m->ILE-1;i<b_a_m->ITE;i++)
//     x[i] = (i - b_a_m->ILE)*delta_x;

//   for(int i=b_a_m->ITE;i<b_a_m->IMAX;i++)
//     x[i] = x[i-1] + (x[i-1] - x[i-2])*b_a_m->XSF;

//   for(int i=b_a_m->ILE-2;i>=0;i--)
//     x[i] = x[i+1] + (x[i+1] - x[i+2])*b_a_m->XSF;

//   y[0] = -delta_x/2.;
//   y[1] = -y[0];

//   for(int j=2;j<b_a_m->JMAX;j++)
//     y[j] = y[j-1] + (y[j-1] - y[j-2])*b_a_m->YSF;
// }

// /*
// - gigiaero, 23/03/2026, 1021 hours
// */
// double bi_air_shape(double t,double x_i){
//   return 2.*t - 4.*t*x_i;
// }

// /*
// - gigiaero, 23/03/2026, 1022 hours
// */
// void bi_air_dirichlet_vals_free(double *x,b_conditions_2d *b_c,double uinf){
//   for(int i=b_c->range[0],idx=0;i<=b_c->range[1];i++,idx++)
//     b_c->val[idx] = uinf*x[i];
// }

// /*
// - gigiaero, 23/03/2026, 1022 hours
// */
// void bi_air_dirichlet_vals_wall(double *x,b_conditions_2d *b_c,double uinf,
//                                 double t){
//   for(int i=b_c->range[0],idx=0;i<=b_c->range[1];i++,idx++)
//     b_c->val[idx] = uinf*bi_air_shape(t,x[i]);
// }




int main(){

  // Solution configurations
  sim_parameters config;
  config.Ntype = 2;
  config.max_iter = 5000;
  config.qtimes = 10;
  config.save_last_only = 1;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/bi_air";

  // Mesh
  bi_air_mesh b_a_mesh;
  b_a_mesh.ILE = 11;   // Leading edge
  b_a_mesh.ITE = 31;   // Trailing edge
  b_a_mesh.IMAX = 41;  // Number of points along x
  b_a_mesh.JMAX = 12;  // Number of points along y
  b_a_mesh.XSF = 1.25; // Stretching factor, x
  b_a_mesh.YSF = 1.25; // Stretching factor, y

  // Physical properties
  double uinf = 1.;
  double t = 0.05;
  // double alpha = 1.; // CONFIRMAR ISTO
  // double lx = 4.;
  // double ly = 2.;

  // Initialize mesh and temperature array
  double *x = malloc(sizeof(double)*b_a_mesh.IMAX);
  double *y = malloc(sizeof(double)*b_a_mesh.JMAX);
  double (*phi)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *phi);
  biconvex_airfoil_mesh(&b_a_mesh,x,y);

  // print_1d_array(b_a_mesh.IMAX,x);
  // putchar('\n');
  // print_1d_array(b_a_mesh.JMAX,y);


  // Boundary conditions
  int num_b_c = 6;
  b_conditions_2d b_c[num_b_c];
  // Down, airfoil
  b_c[0].type = 'd'; 
  b_c[0].axis = 1;
  b_c[0].position = 0;
  b_c[0].range[0] = b_a_mesh.ILE - 1;
  b_c[0].range[1] = b_a_mesh.ITE - 1;
  b_c[0].val = malloc(sizeof(double)*(b_c[0].range[1] - b_c[0].range[0] + 1));
  bi_air_dirichlet_vals_wall(x,&b_c[0],uinf,t);
  // Down, upstream of the leading edge
  b_c[1].type = 'd'; 
  b_c[1].axis = 1;
  b_c[1].position = 0;
  b_c[1].range[0] = 0;
  b_c[1].range[1] = b_a_mesh.ILE - 2;
  b_c[1].val = malloc(sizeof(double)*(b_c[1].range[1] - b_c[1].range[0] + 1));
  bi_air_dirichlet_vals_free(x,&b_c[1],uinf);
  // Down, downstream of the leading edge
  b_c[2].type = 'd'; 
  b_c[2].axis = 1;
  b_c[2].position = 0;
  b_c[2].range[0] = b_a_mesh.ITE;
  b_c[2].range[1] = b_a_mesh.IMAX - 1;
  b_c[2].val = malloc(sizeof(double)*(b_c[2].range[1] - b_c[2].range[0] + 1));
  bi_air_dirichlet_vals_free(x,&b_c[2],uinf);
  // Left (inlet)
  b_c[3].type = 'D';
  b_c[3].axis = 2;
  b_c[3].position = 0;
  b_c[3].range[0] = 1;
  b_c[3].range[1] = b_a_mesh.JMAX - 1;
  b_c[3].val = malloc(sizeof(double));
  b_c[3].val[0] = uinf*x[0];
  // Right (outflow)
  b_c[4].type = 'D'; 
  b_c[4].axis = 2;
  b_c[4].position = b_a_mesh.IMAX - 1;
  b_c[4].range[0] = 1;
  b_c[4].range[1] = b_a_mesh.JMAX - 1;
  b_c[4].val = malloc(sizeof(double));
  b_c[4].val[0] = uinf*x[b_a_mesh.IMAX-1];
  // Up (open boundary)
  b_c[5].type = 'd'; 
  b_c[5].axis = 1;
  b_c[5].position = b_a_mesh.JMAX - 1;
  b_c[5].range[0] = 1;
  b_c[5].range[1] = b_a_mesh.IMAX - 2;
  b_c[5].val = malloc(sizeof(double)*(b_c[5].range[1] - b_c[5].range[0] + 1));
  bi_air_dirichlet_vals_free(x,&b_c[5],uinf);
  
  // Reapplied boundary conditions
  int num_b_c_r = 3;
  b_conditions_2d b_c_r[num_b_c_r];
  // Down, airfoil
  b_c[0].type = 'd'; 
  b_c[0].axis = 1;
  b_c[0].position = 0;
  b_c[0].range[0] = b_a_mesh.ILE - 1;
  b_c[0].range[1] = b_a_mesh.ITE - 1;
  b_c[0].val = malloc(sizeof(double)*(b_c[0].range[1] - b_c[0].range[0] + 1));
  bi_air_dirichlet_vals_wall(x,&b_c[0],uinf,t);
  // Down, upstream of the leading edge
  b_c[1].type = 'd'; 
  b_c[1].axis = 1;
  b_c[1].position = 0;
  b_c[1].range[0] = 0;
  b_c[1].range[1] = b_a_mesh.ILE - 2;
  b_c[1].val = malloc(sizeof(double)*(b_c[1].range[1] - b_c[1].range[0] + 1));
  bi_air_dirichlet_vals_free(x,&b_c[1],uinf);
  // Down, downstream of the leading edge
  b_c[2].type = 'd'; 
  b_c[2].axis = 1;
  b_c[2].position = 0;
  b_c[2].range[0] = b_a_mesh.ITE;
  b_c[2].range[1] = b_a_mesh.IMAX - 1;
  b_c[2].val = malloc(sizeof(double)*(b_c[2].range[1] - b_c[2].range[0] + 1));
  bi_air_dirichlet_vals_free(x,&b_c[2],uinf);



  // Initialize boundary conditions
  apply_b_c(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,num_b_c,b_c);

  // Initialize initial conditions
  for(int j=1;j<b_a_mesh.JMAX-1;j++)
    copy_1d_array_range(1,b_a_mesh.IMAX-1,phi[b_a_mesh.JMAX-1],phi[j]);

  // Solve
  config.casename = malloc(sizeof(char)*200);
  strcpy(config.casename,output_file);
  // evaluate_delta_form(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,x,y,&config);

  save_mesh(b_a_mesh.JMAX,b_a_mesh.IMAX,x,y,config.casename);

  char filename[] = "test.txt";
  print_2d_array_to_file(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,filename,1);




  for(int i=0;i<num_b_c;i++)
    free(b_c[i].val);
  free(x);
  free(y);
  free(phi);
  free(config.casename);

  return 0;
}