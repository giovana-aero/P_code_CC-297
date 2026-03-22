#include<math.h>
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
- gigiaero, 22/03/2026, 0958 hours
*/
void calc_residual(double *phi_elem,double *phi_old_elem,double *N_elem,
                   double *res_elem){

  double val;
  val = fabs((*N_elem)*((*phi_elem) - (*phi_old_elem)));

  if(val > *res_elem)
    *res_elem = val;
}

/*
- gigiaero, 22/03/2026, 1034 hours
*/
void save_mesh(int m,int n,double *x,double *y,char *casename){
  char *filename = malloc(sizeof(char)*200);
  int str_end_idx;

  find_str_end(filename,&str_end_idx);

  strcat(filename,casename);
  find_str_end(filename,&str_end_idx);
  strcat(filename,"_mesh_x.msh");
  print_1d_array_to_file(n,x,filename);
  filename[str_end_idx] = '\0';
  strcat(filename,"_mesh_y.msh");
  print_1d_array_to_file(m,y,filename);

  free(filename);

}

/*
- gigiaero, 22/03/2026, 1121 hours
*/
void save_results_qtimes(int m,int n,double phi[m][n],int *iter,int *qtimes,
                         char *buffer,char *filename_save,int *str_end_idx){
  if(*iter%(*qtimes) == 0){
    sprintf(buffer,"%010d",*iter);
    strcat(filename_save,buffer);
    strcat(filename_save,".dat");
    print_2d_array_to_file(m,n,phi,filename_save,1);
    filename_save[*str_end_idx] = '\0'; // reset string to replace iter number
  }
}

/*
- gigiaero, 19/03/2026, 2338 hours
*/
void solve_p_jacobi_2d_rectangular(int m,int n,double phi[m][n],double *x,
                                   double *y,sim_parameters *config){
  // Solver variables
  double (*phi_old)[n] = calloc(m,sizeof *phi_old);
  double (*N)[n] = calloc(m,sizeof *N);
  double L_phi;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int q_counter = 0, str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double *res = calloc(config->max_iter,sizeof(double));
  // double val;
  FILE *file_log;

  // Configure log file
  strcat(filename_log,config->casename);
  strcat(filename_log,".log");
  // Reset it, if exists, then reopen
  file_log = fopen(filename_log,"w");
  fclose(file_log);
  file_log = fopen(filename_log,"a");

  // Save mesh
  save_mesh(m,n,x,y,config->casename);
  strcat(filename_save,config->casename);

  // Prepare string to save simulation data  
  strcat(filename_save,"_iter_");
  find_str_end(filename_save,&str_end_idx);

  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++)
      N_p_jacobi(&N[j][i],x,y,i,j);
  }

  for(int iter=1;iter<=config->max_iter;iter++){
    copy_2d_array(m,n,phi,phi_old);

    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi,m,n,phi,x,y,i,j);
        phi[j][i] = -L_phi/N[j][i] + phi_old[j][i];
        calc_residual(&phi[j][i],&phi_old[j][i],&N[j][i],&res[iter]);
      }
    }

    printf("Iteration %010d | Res %.6e\n",iter,res[iter]);

    save_results_qtimes(m,n,phi,&iter,&(config->qtimes),buffer,filename_save,&str_end_idx);

    fprintf(file_log,"%.6e\n",res[iter]);

    if(res[iter] <= config->eps){
      puts("<< Convergence! >>");
      break;
    }
  }

  free(phi_old);
  free(N);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
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
  config.max_iter = 50;
  config.qtimes = 10;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/laplace2d";

  // Mesh
  int m = 50;
  int n = 50;

  // Physical properties
  // double alpha = 1.; // CONFIRMAR ISTO
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

  // Initialize boundary conditions
  apply_b_cs(m,n,T,num_b_cs,b_c);

  // Solve
  config.casename = malloc(sizeof(char)*200);
  strcat(config.casename,output_file);
  evaluate_delta_form(m,n,T,x,y,&config);

  free(config.casename);

  return 0;
}