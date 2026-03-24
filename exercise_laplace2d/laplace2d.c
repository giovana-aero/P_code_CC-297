#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/mesh.h"
#include"../include/num_methods.h"
#include"../include/strings.h"



void solve_lgs_rectangular(int m,int n,double phi[m][n],double *x,
                              double *y,sim_parameters *config,
                              int num_b_c_r,b_conditions_2d *b_c_r){
  // Solver variables
  double (*phi_old)[n] = calloc(m,sizeof *phi_old);
  // double (*N)[n-2] = calloc(m-2,sizeof *N);
  double *dx2 = malloc(sizeof(double)*n);
  // double *dy2 = malloc(sizeof(double)*m);
  double L_phi,vars_old,val;
  double r2 = 2./config->r;
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double *res = calloc(config->max_iter,sizeof(double));
  FILE *file_log;

  // Configure log file
  strcpy(filename_log,config->casename);
  strcat(filename_log,".log");
  // Reset it, if exists, then reopen
  file_log = fopen(filename_log,"w");
  fclose(file_log);
  file_log = fopen(filename_log,"a");

  // Prepare string to save simulation data  
  strcpy(filename_save,config->casename);
  strcat(filename_save,"_iter_");
  find_str_end(filename_save,&str_end_idx);

  for(int i=0;i<n;i++){
    dx2[i] = delta_xy(x,i);
    dx2[i] *= dx2[i];
  }

  // for(int j=0;j<m;j++){
  //   dy2[j] = delta_xy(y,j);
  //   dy2[j] *= dy2[j];
  // }

  // Save initial condition
  save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    copy_2d_array(m,n,phi,phi_old);

    // Recalculate boundary conditions (repeated b_cs)
    apply_b_c(m,n,phi,num_b_c_r,b_c_r,x,y); 

    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi,m,n,phi,x,y,i,j);
        vars_old = -L_phi + (phi_old[j][i-1] - r2*phi_old[j][i])/dx2[i] + 
                   (phi_old[j-1][i] - r2*phi_old[j][i])/dy2[j];
        phi[j][i] = (vars_old*dx2[i]*dy2[j] - dy2[j]*phi[j][i-1] - 
                     dx2[i]*phi[j-1][i])/(-r2*(dx2[i] + dy2[j]));

        val = (phi[j][i-1] - r2*phi[j][i])/dx2[i] + 
              (phi[j-1][i] - r2*phi[j][i])/dy2[j] - 
              (phi_old[j][i-1] - r2*phi_old[j][i])/dx2[i] - 
              (phi_old[j-1][i] - r2*phi_old[j][i])/dy2[j];
        val = fabs(val);
        if(val > res[iter])
          res[iter] = val;
      }
    } 

    printf("Iteration %010d | Res %.6e\n",iter,res[iter]);

    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

    fprintf(file_log,"%.6e\n",res[iter]);

    if(res[iter] >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    if(res[iter] <= config->eps){
      puts("<< Convergence! >>");
      iter++;
      break;
    }
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0){
    // buffer[0] = '\0';
    sprintf(buffer,"L");
    config->save_last_only = 0;
    iter--;
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
    &str_end_idx);
  }

  free(phi_old);
  free(dx2);
  // free(dy2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}



int main(){

  // Solution configurations
  sim_parameters config;
  config.Ntype = 3;
  config.r = .95;
  config.max_iter = 50000;
  config.qtimes = 10;
  config.save_last_only = 1;
  config.eps = 1.e-5; // Convergence criterion
  char output_file[] = "results/laplace2d";

  // Mesh
  int m = 100;
  int n = 100;

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
  b_c[counter].val[0] = 10.;
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
  b_c[counter].val[0] = 0.;
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