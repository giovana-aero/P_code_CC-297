#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/bi_air_lib.h"
#include"../include/mesh.h"
#include"../include/num_methods.h"
#include"../include/strings.h"

#define div_ref 1e100

/*
- gigiaero, 22/03/2026, 2131 hours
*/
void biconvex_airfoil_mesh(bi_air_phys_mesh *b_a_m,double *x, double*y){

  double delta_x = 1./(b_a_m->ITE -b_a_m->ILE);

  // Airfoil
  for(int i=b_a_m->ILE-1;i<b_a_m->ITE;i++)
    x[i] = ((i + 1.) - b_a_m->ILE)*delta_x;

  // Downstream
  for(int i=b_a_m->ITE;i<b_a_m->IMAX;i++)
    x[i] = x[i-1] + (x[i-1] - x[i-2])*b_a_m->XSF;

  // Upstream
  for(int i=b_a_m->ILE-2;i>=0;i--)
    x[i] = x[i+1] + (x[i+1] - x[i+2])*b_a_m->XSF;

  y[0] = -delta_x/2.;
  y[1] = -y[0];

  for(int j=2;j<b_a_m->JMAX;j++)
    y[j] = y[j-1] + (y[j-1] - y[j-2])*b_a_m->YSF;
}

/*
- gigiaero, 23/03/2026, 1021 hours
*/
double bi_air_shape(double t,double x_i){
  return 2.*t - 4.*t*x_i;
}

/*
- gigiaero, 23/03/2026, 1022 hours
*/
void bi_air_dirichlet_vals_free(double *x,b_conditions_2d *b_c,double uinf){
  for(int i=b_c->range[0],idx=0;i<=b_c->range[1];i++,idx++)
    b_c->val[idx] = uinf*x[i];
}

/*
- gigiaero, 23/03/2026, 1022 hours
*/
void bi_air_dirichlet_vals_wall(double *x,b_conditions_2d *b_c,double uinf,
                                double t){
  for(int i=b_c->range[0],idx=0;i<=b_c->range[1];i++,idx++)
    b_c->val[idx] = uinf*bi_air_shape(t,x[i]);
}

/*
- gigiaero, 25/03/2026, 1521 hours
*/
void evaluate_delta_form_bi_air(int m,int n,double phi[m][n],double *x,
                                double *y,sim_parameters *config,
                                bi_air_phys_mesh *b_a_m){
  save_mesh(m,n,x,y,config->casename);

  switch(config->Ntype){
    case 1:
      puts("Point Jacobi, biconvex airfoil");
      solve_p_jacobi_2d_rectangular_bi_air(m,n,phi,x,y,config,b_a_m);
      break;

    case 2:
      puts("Gauss-Seidel, biconvex airfoil");
      // solve_g_seidel_2d_rectangular_bi_air(m,n,phi,x,y,config,num_b_c_r,b_c_r);
      break;

    case 3:
      puts("SOR, biconvex airfoil");
      // solve_sor_2d_rectangular_bi_air(m,n,phi,x,y,config,num_b_c_r,b_c_r);
      break;

    case 4:
      puts("Line-Gauss-Seidel, biconvex airfoil");
      // solve_lgs_2d_rectangular_bi_air(m,n,phi,x,y,config,num_b_c_r,b_c_r);
      break;

    case 5:
      puts("SLOR, biconvex airfoil");
      // solve_slor_2d_rectangular_bi_air(m,n,phi,x,y,config,num_b_c_r,b_c_r);
      break;

    default:
      puts("evaluate_delta_form: Invalid Ntype");
      exit(32);
  }
}

/*
- gigiaero, 25/03/2026, 1606 hours
*/
void get_u_v_potential(int m,int n,double phi[m][n],double u[m][n],
                       double v[m][n],double Ve[m][n],double *x,double *y){
  // Domain interior
  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++){
      scheme_der1_o2_central(&u[j][i],m,n,phi,x,i,j,1);
      scheme_der1_o2_central(&v[j][i],m,n,phi,y,i,j,2);
    }
  }

  // Vertical boundaries
  for(int j=0;j<m;j++){
    scheme_der1_o2_forward(&u[j][0],m,n,phi,x,0,j,1);
    scheme_der1_o2_backward(&u[j][n-1],m,n,phi,x,n-1,j,1);
    if(j>0 && j<m-1){
      scheme_der1_o2_central(&v[j][0],m,n,phi,y,0,j,2);
      scheme_der1_o2_central(&v[j][n-1],m,n,phi,y,n-1,j,2);
    }
  }

  // Horizontal boundaries
  for(int i=0;i<n;i++){
    scheme_der1_o2_forward(&v[0][i],m,n,phi,y,i,0,2);
    scheme_der1_o2_backward(&v[m-1][i],m,n,phi,y,i,m-1,2);
    if(i>0 && i<n-1){
      scheme_der1_o2_central(&u[0][i],m,n,phi,x,i,0,1);
      scheme_der1_o2_central(&u[m-1][i],m,n,phi,x,i,m-1,1);
    }
  }

  // Velocity resultant
  for(int j=m-1;j>=0;j--){
    for(int i=0;i<n;i++)
    Ve[j][i] = sqrt(pow(u[j][i],2) + pow(v[j][i],2));
  }
}

/*
- gigiaero, 25/03/2026, 1521 hours
*/
void solve_p_jacobi_2d_rectangular_bi_air(int m,int n,double phi[m][n],
                                          double *x,double *y,
                                          sim_parameters *config,
                                          bi_air_phys_mesh *b_a_m){
  // Solver variables
  double (*phi_old)[n] = calloc(m,sizeof *phi_old);
  double (*N)[n-2] = calloc(m-1,sizeof *N);
  double *tmp_y = malloc(sizeof(double)*3);
  double (*tmp_phi)[3] = calloc(3,sizeof *tmp_phi);
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

  tmp_y[0] = -y[2];
  tmp_y[1] = y[0];
  tmp_y[2] = y[1];
  for(int i=1;i<n-1;i++){
    N_p_jacobi(&N[0][i-1],x,tmp_y,i,1);
    for(int j=1;j<m-1;j++)
      N_p_jacobi(&N[j][i-1],x,y,i,j);
  }

  // Save initial condition
  save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  // print_2d_array(m,n,phi);
  // putchar('\n');
  // char filename[] = "test3.txt";

  for(iter;iter<=config->max_iter;iter++){
    copy_2d_array(m,n,phi,phi_old);

    for(int i=1;i<n-1;i++){
      // printf("<< %d >>\n",i);
      if(i >= b_a_m->ILE-1 && i < b_a_m->ITE){  // Airfoil
        // printf("%02d | airfoil\n",i);
        phi[0][i] = phi[1][i] - (y[1] - y[0])*b_a_m->uinf*bi_air_shape(b_a_m->t,x[i]);
        // printf("%f | %f | %f\n",(y[1] - y[0]),b_a_m->uinf,bi_air_shape(b_a_m->t,x[i]));
      }
      else{  // Free stream
        // printf("%02d | free stream\n",i);
        build_tmp_A_neumann_y_down(m,n,phi_old,tmp_phi,tmp_y,0.,i);
        // tmp_phi[1][0] = phi[0][i-1];
        // tmp_phi[1][1] = phi[0][i];
        // tmp_phi[1][2] = phi[0][i+1];
        // tmp_phi[2][1] = phi[1][i];
        // tmp_phi[0][1] = phi[1][i] + 2.*delta_xy(tmp_y,1)*0.;
        // printf("%f\n",phi[0][i-1]);
        // printf("%f\n",phi[0][i]);
        // printf("%f\n",phi[0][i+1]);
        // printf("%f\n",phi[0][i] + 2.*delta_xy(tmp_y,1)*0.);
        scheme_der2_o2_central_var_deltas_xy(&L_phi,3,3,tmp_phi,x,tmp_y,1,1);
        phi[0][i] = -L_phi/N[0][i-1] + phi[0][i];

        // printf("%f | %f | %f\n",L_phi,N[0][i-1],phi[0][i]);

        // print_2d_array(3,3,tmp_phi);
        // print_1d_array(3,tmp_y);
        // putchar('\n');
      }
    }

    // print_2d_array_to_file(m,n,phi,filename,1);
    // return;

    for(int i=1;i<n-1;i++){
      for(int j=1;j<m-1;j++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi,m,n,phi_old,x,y,i,j);
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
    sprintf(buffer,"L");
    config->save_last_only = 0;
    iter--;
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                        &str_end_idx);
  }

  free(phi_old);
  free(N);
  free(tmp_phi);
  free(tmp_y);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}