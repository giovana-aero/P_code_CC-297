#include<math.h>
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
  config.r = 1.95;
  config.max_iter = 10000;
  config.qtimes = 1;
  config.save_i_c = 0;
  config.save_last_only = 1;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/bi_air";

  // Mesh
  bi_air_phys_mesh b_a_mesh;
  // // ----------------------------------------------- original configuration
  b_a_mesh.ILE = 11;   // Leading edge
  b_a_mesh.ITE = 31;   // Trailing edge
  b_a_mesh.IMAX = 41;  // Number of points along x
  b_a_mesh.JMAX = 12;  // Number of points along y
  b_a_mesh.XSF = 1.25; // Stretching factor, x
  b_a_mesh.YSF = 1.25; // Stretching factor, y
  // ----------------------------------------------- 2x
  // b_a_mesh.ILE = 21;   // Leading edge
  // b_a_mesh.ITE = 61;   // Trailing edge
  // b_a_mesh.IMAX = 82;  // Number of points along x
  // b_a_mesh.JMAX = 24;  // Number of points along y
  // b_a_mesh.XSF = 1.05; // Stretching factor, x
  // b_a_mesh.YSF = 1.05; // Stretching factor, y
  // ----------------------------------------------- 
  // b_a_mesh.ILE = ;   // Leading edge
  // b_a_mesh.ITE = ;   // Trailing edge
  // b_a_mesh.IMAX = ;  // Number of points along x
  // b_a_mesh.JMAX = ;  // Number of points along y
  // b_a_mesh.XSF = 1.25; // Stretching factor, x
  // b_a_mesh.YSF = 1.25; // Stretching factor, y
  // ----------------------------------------------- good example for p_jacobi (?)
  // b_a_mesh.ILE = 21;   // Leading edge
  // b_a_mesh.ITE = 61;   // Trailing edge
  // b_a_mesh.IMAX = 82;  // Number of points along x
  // b_a_mesh.JMAX = 24;  // Number of points along y
  // b_a_mesh.XSF = 1.05; // Stretching factor, x
  // b_a_mesh.YSF = 1.05; // Stretching factor, y

  int chord_l = b_a_mesh.ITE - b_a_mesh.ILE + 1;

  // Physical properties
  b_a_mesh.t = 0.05;
  b_a_mesh.uinf = 1.;

  // Initialize mesh and velocity potential array
  double *x = malloc(sizeof(double)*b_a_mesh.IMAX);
  double *y = malloc(sizeof(double)*b_a_mesh.JMAX);
  double (*phi)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *phi);
  double (*u)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *u);
  double (*v)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *v);
  double (*Ve)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *Ve);
  double (*cp)[b_a_mesh.IMAX] = calloc(b_a_mesh.JMAX,sizeof *cp);
  double *cp_chord = malloc(sizeof(double)*chord_l);
  biconvex_airfoil_mesh(&b_a_mesh,x,y);
  // uniform_rectangular_mesh(b_a_mesh.JMAX,b_a_mesh.IMAX,100,10,x,y);

  // Boundary conditions
  int num_b_c = 3;
  b_conditions_2d b_c[num_b_c];
  int counter = 0;
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
  counter++;
    // // Up (test)
    // b_c[counter].type = 'D'; 
    // b_c[counter].axis = 1;
    // b_c[counter].position = b_a_mesh.JMAX - 1;
    // b_c[counter].range[0] = 1;
    // b_c[counter].range[1] = b_a_mesh.IMAX - 2;
    // b_c[counter].val = malloc(sizeof(double));
    // b_c[counter].val[0] = 0.;
    // counter++;
  // // Down (test)
  // b_c[counter].type = 'D'; 
  // b_c[counter].axis = 1;
  // b_c[counter].position = b_a_mesh.JMAX - 1;
  // b_c[counter].range[0] = 1;
  // b_c[counter].range[1] = b_a_mesh.IMAX - 2;
  // b_c[counter].val = malloc(sizeof(double));
  // b_c[counter].val[0] = 0.;

  // Reapplied boundary conditions
  // int num_b_c_r = 0;

  // Initialize boundary conditions and initial conditions
  apply_b_c(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,num_b_c,b_c,x,y);
  for(int j=1;j<b_a_mesh.JMAX-1;j++)
    copy_1d_array_range(1,b_a_mesh.IMAX-1,phi[b_a_mesh.JMAX-1],phi[j]);

  apply_down_b_c(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,x,y,&b_a_mesh);

  // Solve
  config.casename = malloc(sizeof(char)*200);
  strcpy(config.casename,output_file);
  evaluate_delta_form_bi_air(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,x,y,&config,
                             &b_a_mesh);
  // evaluate_delta_form(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,x,y,&config,num_b_c_r,NULL);

  // Post-processing: get velocities and coefficients of pressure
  get_u_v_potential(b_a_mesh.JMAX,b_a_mesh.IMAX,phi,u,v,Ve,x,y);
  get_cp_bi_air(b_a_mesh.JMAX,b_a_mesh.IMAX,cp,u,v,x,y,&b_a_mesh);
  get_cp_bi_air_chord(cp_chord,b_a_mesh.JMAX,b_a_mesh.IMAX,phi,u,x,y,&b_a_mesh);

  // Save everything
  char *save_files = malloc(sizeof(char)*200);
  sprintf(save_files,"%s_u.dat",config.casename);
  print_2d_array_to_file(b_a_mesh.JMAX,b_a_mesh.IMAX,u,save_files,0);
  sprintf(save_files,"%s_v.dat",config.casename);
  print_2d_array_to_file(b_a_mesh.JMAX,b_a_mesh.IMAX,v,save_files,0);
  sprintf(save_files,"%s_Ve.dat",config.casename);
  print_2d_array_to_file(b_a_mesh.JMAX,b_a_mesh.IMAX,Ve,save_files,0);
  sprintf(save_files,"%s_cp.dat",config.casename);
  print_2d_array_to_file(b_a_mesh.JMAX,b_a_mesh.IMAX,cp,save_files,0);
  sprintf(save_files,"%s_cp_chord.dat",config.casename);
  print_1d_array_to_file(chord_l,cp_chord,save_files);


  for(int i=0;i<num_b_c;i++)
    free(b_c[i].val);
  free(x);
  free(y);
  free(phi);
  free(u);
  free(v);
  free(Ve);
  free(cp);
  free(cp_chord);
  free(config.casename);
  free(save_files);

  return 0;
}