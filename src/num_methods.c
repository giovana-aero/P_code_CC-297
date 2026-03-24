#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/mesh.h"
#include"../include/num_methods.h"
#include"../include/strings.h"

#define div_ref 1e100


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
- gigiaero, 22/03/2026, 1063 hours
*/
int check_num_digits_int(int *num){
  int i = 1;
  double div = 10.;
  double d_num = (double) *num;

  while(d_num/div >= 1){
    div *= 10.;
    i++;
  }

  return i;
}

/*
- gigiaero, 19/03/2026, 2318 hours
*/
double delta_xy(double *xy,int i){
  return (xy[i+1] - xy[i-1])/2.;
}

/*
This was based on series 02's question 5. Also works if some of the diagonals 
are zero (except the main one)
- gigiaero, 18/03/2026, 1620 hours
*/
void diagonal_matrix_solver(int n,double A[n][n],double *f,double *u){

  // Original matrix
  double *a = (double *)malloc(sizeof(double)*(n-1));
  double *b = (double *)malloc(sizeof(double)*n);
  double *c = (double *)malloc(sizeof(double)*(n-1));
  double *d = (double *)malloc(sizeof(double)*(n-2));
  // Decomposed matrices
  double *k = (double *)malloc(sizeof(double)*(n-1));
  double *w = (double *)malloc(sizeof(double)*n);
  double *x = (double *)malloc(sizeof(double)*(n-1));
  double *z = (double *)malloc(sizeof(double)*(n-2));
  // Intermediary solution
  double *y = (double *)malloc(sizeof(double)*n);
  

  if(a == NULL){
    puts("diagonal_matrix_solver: Allocation failed (a)");
    exit(1);
  }
  if(b == NULL){
    puts("diagonal_matrix_solver: Allocation failed (b)");
    exit(1);
  }
  if(c == NULL){
    puts("diagonal_matrix_solver: Allocation failed (c)");
    exit(1);
  }
  if(d == NULL){
    puts("diagonal_matrix_solver: Allocation failed (d)");
    exit(1);
  }
  if(k == NULL){
    puts("diagonal_matrix_solver: Allocation failed (k)");
    exit(1);
  }
  if(w == NULL){
    puts("diagonal_matrix_solver: Allocation failed (w)");
    exit(1);
  }
  if(x == NULL){
    puts("diagonal_matrix_solver: Allocation failed (x)");
    exit(1);
  }
  if(z == NULL){
    puts("diagonal_matrix_solver: Allocation failed (z)");
    exit(1);
  }
  if(y == NULL){
    puts("diagonal_matrix_solver: Allocation failed (y)");
    exit(1);
  }

  // Get coefficients of the original matrix
  for(int i=0;i<n;i++){
    b[i] = A[i][i];

    if(i < n-1)
      a[i] = A[i+1][i];

    if(i > 0)
      c[i-1] = A[i-1][i];

    if(i > 1)
      d[i-2] = A[i-2][i];
  }

  // Calculate coefficients of decomposed matrices
  k[0] = a[0];
  w[0] = b[0];
  x[0] = c[0]/w[0];
  z[0] = d[0]/w[0];
  for(int i=1;i<n;i++){
    if(i < (n - 1))
      k[i] = a[i];

    w[i] = b[i] - k[i-1]*x[i-1];

    if(i < (n - 1))
      x[i] = (c[i] - k[i-1]*z[i-1])/w[i];
    
    if(i < (n - 2))
      z[i] = d[i]/w[i];
  }

  // Solve for y
  y[0] = f[0]/w[0];
  for(int i=1;i<n;i++)
    y[i] = (f[i] - k[i-1]*y[i-1])/w[i];

  // Solve for u
  u[n-1] = y[n-1];
  u[n-2] = y[n-2] - x[n-2]*u[n-1];
  for(int i=n-3;i>=0;i--)
    u[i] = y[i] - x[i]*u[i+1] - z[i]*u[i+2];

  free(a);
  free(b);
  free(c);
  free(d);
  free(k);
  free(w);
  free(x);
  free(z);
  free(y);
}

/*
- gigiaero, 19/03/2026, 2318 hours
*/
void evaluate_delta_form(int m,int n,double phi[m][n],double *x,double *y,
                        sim_parameters *config,int num_b_c_r,
                        b_conditions_2d *b_c_r){
  save_mesh(m,n,x,y,config->casename);

  switch(config->Ntype){
    case 1:
      puts("Point Jacobi");
      solve_p_jacobi_2d_rectangular(m,n,phi,x,y,config,num_b_c_r,b_c_r);
      break;

    case 2:
      puts("Gauss-Seidel");
      solve_g_seidel_2d_rectangular(m,n,phi,x,y,config,num_b_c_r,b_c_r);
      break;

    case 3:
      puts("SOR");
      solve_sor_2d_rectangular(m,n,phi,x,y,config,num_b_c_r,b_c_r);
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

/*
- gigiaero, 19/03/2026, 2318 hours
*/
void N_p_jacobi(double *N,double *x,double *y,int i,int j){
  *N = -(2./(delta_xy(x,i)*delta_xy(x,i)) + 
         2./(delta_xy(y,j)*delta_xy(y,j)));
}

/*
- gigiaero, 22/03/2026, 1121 hours

in this order, the conditions on the internal "if" represent the normal 
operation of the function, the first iteration, and the last iteration
- gigiaero, 22/03/2026 hours
*/
void save_results_qtimes(int m,int n,double phi[m][n],int *iter,
                        sim_parameters *s_p,char *buffer,char *filename_save,
                        int *str_end_idx){
  if(!s_p->save_last_only){
    if((*iter)%(s_p->qtimes) == 0 || (*iter) == 0 || buffer[0] == 'L'){
      sprintf(buffer,"%010d",*iter);
      strcat(filename_save,buffer);
      strcat(filename_save,".dat");
      print_2d_array_to_file(m,n,phi,filename_save,0);
      filename_save[*str_end_idx] = '\0'; // reset string to replace iter number

      if((*iter) == 0)
        (*iter)++;
    }
  }
}


/*
Scheme: second derivative, second order, central
- gigiaero, 13/03/2026, 2315 hours
*/
double scheme_der2_o2_central(double phi_ip1,double phi_i,double phi_im1,
                              double x_ip1,double x_i,double x_im1){
  return 2./(x_ip1 - x_im1)*((phi_ip1 - phi_i)/(x_ip1 - x_i) - 
         (phi_i - phi_im1)/(x_i - x_im1));
}

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
- gigiaero, 23/03/2026 hours
*/
void solve_g_seidel_2d_rectangular(int m,int n,double phi[m][n],double *x,
                                   double *y,sim_parameters *config,
                                   int num_b_c_r,b_conditions_2d *b_c_r){
  // Solver variables
  double (*phi_old)[n] = calloc(m,sizeof *phi_old);
  // double (*N)[n-2] = calloc(m-2,sizeof *N);
  double *dx2 = malloc(sizeof(double)*n);
  double *dy2 = malloc(sizeof(double)*m);
  double L_phi,vars_old,val;
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

  for(int j=0;j<m;j++){
    dy2[j] = delta_xy(y,j);
    dy2[j] *= dy2[j];
  }

  // Save initial condition
  save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    copy_2d_array(m,n,phi,phi_old);

    // Recalculate boundary conditions (repeated b_cs)
    apply_b_c(m,n,phi,num_b_c_r,b_c_r,x,y); 

    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi,m,n,phi,x,y,i,j);
        vars_old = -L_phi + (phi_old[j][i-1] - 2.*phi_old[j][i])/dx2[i] + 
                   (phi_old[j-1][i] - 2.*phi_old[j][i])/dy2[j];
        phi[j][i] = (vars_old*dx2[i]*dy2[j] - dy2[j]*phi[j][i-1] - 
                     dx2[i]*phi[j-1][i])/(-2.*(dx2[i] + dy2[j]));

        val = (phi[j][i-1] - 2.*phi[j][i])/dx2[i] + 
              (phi[j-1][i] - 2.*phi[j][i])/dy2[j] - 
              (phi_old[j][i-1] - 2.*phi_old[j][i])/dx2[i] - 
              (phi_old[j-1][i] - 2.*phi_old[j][i])/dy2[j];
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
  free(dy2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}

/*
- gigiaero, 19/03/2026, 2338 hours
*/
void solve_p_jacobi_2d_rectangular(int m,int n,double phi[m][n],double *x,
                                   double *y,sim_parameters *config,
                                   int num_b_c_r,b_conditions_2d *b_c_r){
  // Solver variables
  double (*phi_old)[n] = calloc(m,sizeof *phi_old);
  double (*N)[n-2] = calloc(m-2,sizeof *N);
  double L_phi;
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

  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++)
      N_p_jacobi(&N[j-1][i-1],x,y,i,j);
  }

  // Save initial condition
  save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    copy_2d_array(m,n,phi,phi_old);

    // Recalculate boundary conditions (repeated b_cs)
    apply_b_c(m,n,phi,num_b_c_r,b_c_r,x,y); 

    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi,m,n,phi,x,y,i,j);
        phi[j][i] = -L_phi/N[j-1][i-1] + phi_old[j][i];
        calc_residual(&phi[j][i],&phi_old[j][i],&N[j-1][i-1],&res[iter]);
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
  free(N);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}

void solve_sor_2d_rectangular(int m,int n,double phi[m][n],double *x,
                              double *y,sim_parameters *config,
                              int num_b_c_r,b_conditions_2d *b_c_r){
  // Solver variables
  double (*phi_old)[n] = calloc(m,sizeof *phi_old);
  // double (*N)[n-2] = calloc(m-2,sizeof *N);
  double *dx2 = malloc(sizeof(double)*n);
  double *dy2 = malloc(sizeof(double)*m);
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

  for(int j=0;j<m;j++){
    dy2[j] = delta_xy(y,j);
    dy2[j] *= dy2[j];
  }

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
  free(dy2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}