#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/mesh.h"
#include"../include/num_methods.h"
#include"../include/strings.h"

/*
- gigiaero, 19/03/2026, 2333 hours
*/
void scheme_der2_o2_central_var_deltas_xy(double *f,int m,int n,
                                          double phi[m][n],double *x,double *y,
                                          int i,int j){
  *f = 2./(x[i+1] - x[i-1])*
       ((phi[j][i+1] - phi[j][i])/(x[i+1] - x[i]) - 
       (phi[j][i] - phi[j][i-1])/(x[i] - x[i-1])) + 
       2./(y[j+1] - y[j-1])*
       ((phi[j+1][i] - phi[j][i])/(y[j+1] - y[j]) - 
       (phi[j][i] - phi[j-1][i])/(y[j] - y[j-1]));
}

/*
- gigiaero, 19/03/2026, 2318 hours
*/
void N_p_jacobi(double *N,double *x,double *y,int i,int j){
  *N = -(2./(delta_xy(x,i)*delta_xy(x,i)) + 
         2./(delta_xy(y,j)*delta_xy(y,j)));
}

/*
- gigiaero, 19/03/2026, 2338 hours
*/
void solve_p_jacobi_2d_rectangular(int m,int n,double phi[m][n],double *x,
                                   double *y,sim_parameters *config){

  double phi_old[m][n];
  int res[config->max_iter];
  double N;
  char filename_save[200], filename_log[200], buffer[100];
  int q_counter = 0, str_end_idx;

  // Configure log file
  strcat(filename_log,"_log.log");

  // Save mesh
  strcat(filename_save,config->casename);
  find_str_end(filename_save,&str_end_idx);
  strcat(filename_save,"_mesh_x.msh");
  print_1d_array_to_file(n,x,filename_save);
  filename_save[str_end_idx] = '\0';
  strcat(filename_save,"_mesh_y.msh");
  print_1d_array_to_file(m,y,filename_save);

  // Prepare string to save simulation data  
  filename_save[str_end_idx] = '\0';
  strcat(filename_save,"_iter");
  find_str_end(filename_save,&str_end_idx);

  for(int iter=1;iter<=config->max_iter;iter++){
    copy_2d_array(m,n,phi,phi_old);
    double L_phi, N;

    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi,m,n,phi,x,y,i,j);
        N_p_jacobi(&N,x,y,i,j);
        phi[j][i] = -L_phi/N + phi_old[j][i];
      }
    }

    printf("Iteration %d\n",iter);

    if(iter%config->qtimes == 0){
      itoa(iter,buffer,10);
      strcat(filename_save,buffer);
      strcat(filename_save,".dat");
      print_2d_array_to_file(m,n,phi,filename_save,1);
      filename_save[str_end_idx] = '\0'; // reset string to replace iter
    }
  }
}

/*
- gigiaero, 19/03/2026, 2318 hours
*/
void evaluate_delta_form(int m,int n,double phi[m][n],double *x,double *y,
                         sim_parameters *config){

  switch(config->Ntype){
    case 1:
      puts("Point Jacobi");
      solve_p_jacobi_2d_rectangular(m,n,phi,x,y,config);
      break;

    case 2:
      puts("Gauss-Seidel");

      break;

    case 3:
      puts("SOR");

      break;

    case 4:
      puts("Line-Gauss-Seidel");

      break;

    case 5:
      puts("SLOR");

      break;

    default:
      puts("evaluate_delta_form: Invalid Ntype");
      exit(32);
  }
}

int main(){

  // Solution configurations
  sim_parameters config;
  config.Ntype = 1;
  config.max_iter = 100;
  config.qtimes = 10;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/laplace2d";

  // Mesh
  int m = 5;
  int n = 9;

  // Physical properties
  double alpha = 1.; // CONFIRMAR ISTO
  double lx = 4.;
  double ly = 2.;

  // Boundary conditions
  int num_b_cs = 4;
  b_conditions_2d b_c[num_b_cs];
  b_c[0].type = 'D';
  b_c[0].val = 1.;
  b_c[0].axis = 1;
  b_c[0].position = 0;
  b_c[0].range[0] = 1;
  b_c[0].range[1] = n-1;
  b_c[1].type = 'D';
  b_c[1].val = 2.;
  b_c[1].axis = 2;
  b_c[1].position = n-1;
  b_c[1].range[0] = 1;
  b_c[1].range[1] = m-1;
  b_c[2].type = 'D';
  b_c[2].val = 3.;
  b_c[2].axis = 1;
  b_c[2].position = m-1;
  b_c[2].range[0] = 1;
  b_c[2].range[1] = n-1;
  b_c[3].type = 'D';
  b_c[3].val = 4.;
  b_c[3].axis = 2;
  b_c[3].position = 0;
  b_c[3].range[0] = 1;
  b_c[3].range[1] = m-1;

  // Initialize mesh and temperature array
  double x[n], y[m]; 
  double T[m][n];
  double delta_x = lx/(n - 1.);
  double delta_y = ly/(m - 1.);
  zeros_2d_array(m,n,T);
  uniform_rectangular_mesh(m,n,delta_x,delta_y,x,y);
  // print_1d_array(n,x);
  // print_1d_array(m,y);

  // Initialize boundary conditions
  apply_b_cs(m,n,T,num_b_cs,b_c);
  // dirichlet_rectangular_constant(m,n,T,Tu,range_h,m-1,0);
  // dirichlet_rectangular_constant(m,n,T,Td,range_h,0,0);
  // dirichlet_rectangular_constant(m,n,T,Tl,range_v,0,1);
  // dirichlet_rectangular_constant(m,n,T,Tr,range_v,n-1,1);

  // Solve
  strcat(config.casename,output_file);
  // print_2d_array(m,n,T);
  // evaluate_delta_form(m,n,T,x,y,&config);
  // putchar('\n');
  print_2d_array(m,n,T);

  // Print results to file
  // strcat(config.casename,"_test.txt");
  // print_2d_array_to_file(m,n,T,config.casename);

  return 0;
}