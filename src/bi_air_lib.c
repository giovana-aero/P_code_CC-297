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
- gigiaero, 06/04/2026, 1600 hours
*/
void apply_down_b_c(int m,int n,double phi[m][n],double *x,double *y,
                    bi_air_phys_mesh *b_a_m){
  for(int i=1;i<n-1;i++){
    if(i >= b_a_m->ILE-1 && i < b_a_m->ITE){  // Airfoil
      phi[0][i] = phi[1][i] - (y[1] - y[0])*b_a_m->uinf*bi_air_shape_dx(b_a_m->t,x[i]);
    }
    else{  // Free stream
      phi[0][i] = phi[1][i];
    }
  }
}

/*
- gigiaero, 22/03/2026, 2131 hours
*/
void biconvex_airfoil_mesh(bi_air_phys_mesh *b_a_m,double *x, double*y){

  double delta_x = 1./(b_a_m->ITE - b_a_m->ILE);

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

// /*
// - gigiaero, 06/04/2026, 1051 hours
// */
// double bi_air_shape(double t,double x_i){
//   return 2.*t*x_i*(1. - x_i);
// }

/*
- gigiaero, 23/03/2026, 1021 hours
*/
double bi_air_shape_dx(double t,double x_i){
  // return 2.*t*x_i*(1. - x_i);
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
    b_c->val[idx] = uinf*bi_air_shape_dx(t,x[i]);
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
      solve_g_seidel_2d_rectangular_bi_air(m,n,phi,x,y,config,b_a_m);
      break;

    case 3:
      puts("SOR, biconvex airfoil");
      solve_sor_2d_rectangular_bi_air(m,n,phi,x,y,config,b_a_m);
      break;

    case 4:
      puts("Line-Gauss-Seidel, biconvex airfoil");
      solve_lgs_2d_rectangular_bi_air(m,n,phi,x,y,config,b_a_m);
      break;

    case 5:
      puts("SLOR, biconvex airfoil");
      solve_slor_2d_rectangular_bi_air(m,n,phi,x,y,config,b_a_m);
      break;

    default:
      puts("evaluate_delta_form: Invalid Ntype");
      exit(32);
  }
}

/*
- gigiaero, 06/04/2026, 1634 hours
*/
void get_cp_bi_air(int m,int n,double cp[m][n],double u[m][n],double v[m][n],
                   double *x,double *y,bi_air_phys_mesh *b_a_m){
  double uinf2_m1 = 1./(b_a_m->uinf*b_a_m->uinf);
  double u_val2,v_val2;

  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++){
      v_val2 = v[j][i];
      u_val2 = u[j][i];
      v_val2 *= v_val2;
      u_val2 *= u_val2;
      
      cp[j][i] = 1. - (u_val2 + v_val2)*uinf2_m1;
    }
  }
}

/*
- gigiaero, 06/04/2026, 1607 hours
(kudos to gibson helping me find a certain dumb mistake)
*/
void get_cp_bi_air_chord(double *cp,int m,int n,double phi[m][n],double u[m][n],
                         double *x,double *y,bi_air_phys_mesh *b_a_m){
  double uinf2_m1 = 1./(b_a_m->uinf*b_a_m->uinf);
  double u_val2,v_val2;

  for(int i=b_a_m->ILE-1,k=0;i<b_a_m->ITE;i++,k++){
    v_val2 = (phi[1][i] - phi[0][i])/(y[1] - y[0]);
    u_val2 = (u[0][i] + u[1][i])/2.;
    v_val2 *= v_val2;
    u_val2 *= u_val2;
    
    cp[k] = 1. - (u_val2 + v_val2)*uinf2_m1;
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
- gigiaero, 06/04/2026, 1720 hours
*/
void solve_g_seidel_2d_rectangular_bi_air(int m,int n,double phi[m][n],
                                          double *x,double *y,
                                          sim_parameters *config,
                                          bi_air_phys_mesh *b_a_m){
  // Solver variables
  double (*L_phi)[n-2] = calloc(m-2,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double *dx2 = malloc(sizeof(double)*(n-2));
  double *dy2 = malloc(sizeof(double)*(m-2));
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double *res = calloc(config->max_iter + 1,sizeof(double));
  FILE *file_log;

  // Configure log file
  sprintf(filename_log,"%s.log",config->casename);
  file_log = fopen(filename_log,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save,"%s_iter_",config->casename);
  find_str_end(filename_save,&str_end_idx);

  for(int i=1;i<n-1;i++){
    dx2[i-1] = delta_xy(x,i);
    dx2[i-1] *= dx2[i-1];
  }

  for(int j=1;j<m-1;j++){
    dy2[j-1] = delta_xy(y,j);
    dy2[j-1] *= dy2[j-1];
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res[iter])
          res[iter] = fabs(L_phi[j-1][i-1]);
      }
    }
    
    printf("GS Iteration %010d | Res %.6e\n",iter,res[iter]);

    fprintf(file_log,"%.6e\n",res[iter]);

    // Test for convergence
    if(res[iter] <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res[iter] >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        Cij[j][i] = (-L_phi[j-1][i-1] - Cij[j][i-1]/dx2[i-1]
                     - Cij[j-1][i]/dy2[j-1])/
                     (-2./dx2[i-1] - 2./dy2[j-1]);
    }
    
    // Calculate phi^{n+1}
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        phi[j][i] += Cij[j][i];
    }
    
    // Apply boundary conditions (low edge)
    apply_down_b_c(m,n,phi,x,y,b_a_m);

    if(!config->save_last_only)
      save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                          &str_end_idx);
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0 || config->save_last_only){
    iter--; // To get the correct iteration number
    sprintf(buffer,"L");
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                        &str_end_idx);
  }

  free(L_phi);
  free(Cij);
  free(dx2);
  free(dy2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}

/*
- gigiaero , 08/04/2026, 2137 hours
*/
void solve_lgs_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                     double *y,sim_parameters *config,
                                     bi_air_phys_mesh *b_a_m){
  // Solver variables
  double (*L_phi)[n-2] = calloc(m-2,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double (*A)[m-2] = calloc(m-2,sizeof *A);
  double *f = malloc(sizeof(double)*(m-2));
  double *u = malloc(sizeof(double)*(m-2));
  double *dx2 = malloc(sizeof(double)*(n-2));
  double *py1 = malloc(sizeof(double)*(m-2));
  double *py2 = malloc(sizeof(double)*(m-2));
  double vars_old,val;
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double *res = calloc(config->max_iter + 1,sizeof(double));
  FILE *file_log;

  // Configure log file
  sprintf(filename_log,"%s.log",config->casename);
  puts(filename_log);
  file_log = fopen(filename_log,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save,"%s_iter_",config->casename);
  find_str_end(filename_save,&str_end_idx);

  for(int i=1;i<n-1;i++){
    dx2[i-1] = delta_xy(x,i);
    dx2[i-1] *= dx2[i-1];
  }

  for(int j=1;j<m-1;j++){
    py1[j-1] = (y[j+1] - y[j-1])*(y[j+1] - y[j]);
    py2[j-1] = (y[j+1] - y[j-1])*(y[j] - y[j-1]);
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res[iter])
          res[iter] = fabs(L_phi[j-1][i-1]);
      }
    }

    printf("LGS Iteration %010d | Res %.6e\n",iter,res[iter]);

    fprintf(file_log,"%.6e\n",res[iter]);

    // Test for convergence
    if(res[iter] <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res[iter] >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int i=1;i<n-1;i++){
      A[0][0] = -1./dx2[i-1] - 1./py1[0] - 1./py2[0];
      A[0][1] = 1./py1[0];
      f[0] = (-L_phi[0][i-1] - Cij[1][i]/dx2[i-1])/2. - Cij[0][i]/py2[0];

      for(int j=1;j<m-3;j++){
        A[j][j-1] = 1./py2[j];
        A[j][j] = (-1./dx2[i-1] - 1./py1[j] - 1./py2[j]);
        A[j][j+1] = 1./py1[j]/1;
        f[j] = (-L_phi[j][i-1] - Cij[j+1][i]/dx2[i-1])/2.;
      }

      A[m-3][m-4] = 1./py2[m-3];
      A[m-3][m-3] = -1./dx2[i-1] - 1./py1[m-3] - 1./py2[m-3];
      f[m-3] = (-L_phi[m-3][i-1] - Cij[m-2][i]/dx2[i-1])/2. 
               - Cij[m-1][i]/py1[m-3];
      
      diagonal_matrix_solver(m-2,A,f,u);

      for(int j=1;j<m-1;j++)
        Cij[j][i] = u[j-1];
    }

    // Calculate phi^{n+1}
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        phi[j][i] += Cij[j][i];
    }

    // Apply boundary conditions (low edge)
    apply_down_b_c(m,n,phi,x,y,b_a_m);

    if(!config->save_last_only)
      save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                          &str_end_idx);
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0 || config->save_last_only){
    iter--; // To get the correct iteration number
    sprintf(buffer,"L");
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                        &str_end_idx);
  }

  free(L_phi);
  free(Cij);
  free(A);
  free(f);
  free(u);
  free(dx2);
  free(py1);
  free(py2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}

/*
- gigiaero, 06/04/2026, 1602 hours
*/
void solve_p_jacobi_2d_rectangular_bi_air(int m,int n,double phi[m][n],
                                          double *x,double *y,
                                          sim_parameters *config,
                                          bi_air_phys_mesh *b_a_m){
  // Solver variables
  double (*N)[n-2] = calloc(m,sizeof *N-2);
  double (*L_phi)[n] = calloc(m,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  // double res_val = 0.;
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double *res = calloc(config->max_iter + 1,sizeof(double));
  FILE *file_log;

  // Configure log file
  sprintf(filename_log,"%s.log",config->casename);
  file_log = fopen(filename_log,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save,"%s_iter_",config->casename);
  find_str_end(filename_save,&str_end_idx);

  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++)
      N_p_jacobi(&N[j-1][i-1],x,y,i,j);
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j][i],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j][i]) > res[iter])
          res[iter] = L_phi[j][i];
      }
    }
    
    printf("PJ Iteration %010d | Res %.6e\n",iter,res[iter]);

    fprintf(file_log,"%.6e\n",res[iter]);

    // Test for convergence
    if(res[iter] <= config->eps & iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res[iter] >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        Cij[j][i] = -L_phi[j][i]/N[j-1][i-1];
    }
    
    // Calculate phi^{n+1}
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        phi[j][i] += Cij[j][i];
    }
    
    // Apply boundary conditions (low edge)
    apply_down_b_c(m,n,phi,x,y,b_a_m);

    if(!config->save_last_only)
      save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                          &str_end_idx);
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0 || config->save_last_only){
    iter--; // To get the correct iteration number
    sprintf(buffer,"L");
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                        &str_end_idx);
  }

  free(N);
  free(L_phi);
  free(Cij);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}

/*
- gigiaero, 07/04/2026, 2251 hours
*/
void solve_slor_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                      double *y,sim_parameters *config,
                                      bi_air_phys_mesh *b_a_m){
  // Solver variables
  double (*L_phi)[n-2] = calloc(m-2,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double (*A)[m-2] = calloc(m-2,sizeof *A);
  double *f = malloc(sizeof(double)*(m-2));
  double *u = malloc(sizeof(double)*(m-2));
  double *dx2 = malloc(sizeof(double)*(n-2));
  double *py1 = malloc(sizeof(double)*(m-2));
  double *py2 = malloc(sizeof(double)*(m-2));
  double vars_old,val;
  double r = config->r;
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double *res = calloc(config->max_iter + 1,sizeof(double));
  FILE *file_log;

  // Configure log file
  sprintf(filename_log,"%s.log",config->casename);
  puts(filename_log);
  file_log = fopen(filename_log,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save,"%s_iter_",config->casename);
  find_str_end(filename_save,&str_end_idx);

  for(int i=1;i<n-1;i++){
    dx2[i-1] = delta_xy(x,i);
    dx2[i-1] *= dx2[i-1];
  }

  for(int j=1;j<m-1;j++){
    py1[j-1] = (y[j+1] - y[j-1])*(y[j+1] - y[j]);
    py2[j-1] = (y[j+1] - y[j-1])*(y[j] - y[j-1]);
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res[iter])
          res[iter] = fabs(L_phi[j-1][i-1]);
      }
    }

    printf("SLOR Iteration %010d | Res %.6e\n",iter,res[iter]);

    fprintf(file_log,"%.6e\n",res[iter]);

    // Test for convergence
    if(res[iter] <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res[iter] >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int i=1;i<n-1;i++){
      A[0][0] = -1./dx2[i-1] - 1./py1[0] - 1./py2[0];
      A[0][1] = 1./py1[0];
      f[0] = (-L_phi[0][i-1] - Cij[1][i]/dx2[i-1])*r/2. - Cij[0][i]/py2[0];

      for(int j=1;j<m-3;j++){
        A[j][j-1] = 1./py2[j];
        A[j][j] = (-1./dx2[i-1] - 1./py1[j] - 1./py2[j]);
        A[j][j+1] = 1./py1[j]/1;
        f[j] = (-L_phi[j][i-1] - Cij[j+1][i]/dx2[i-1])*r/2.;
      }

      A[m-3][m-4] = 1./py2[m-3];
      A[m-3][m-3] = -1./dx2[i-1] - 1./py1[m-3] - 1./py2[m-3];
      f[m-3] = (-L_phi[m-3][i-1] - Cij[m-2][i]/dx2[i-1])*r/2. 
               - Cij[m-1][i]/py1[m-3];
      
      diagonal_matrix_solver(m-2,A,f,u);

      for(int j=1;j<m-1;j++)
        Cij[j][i] = u[j-1];
    }

    // Calculate phi^{n+1}
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        phi[j][i] += Cij[j][i];
    }

    // Apply boundary conditions (low edge)
    apply_down_b_c(m,n,phi,x,y,b_a_m);

    if(!config->save_last_only)
      save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                          &str_end_idx);
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0 || config->save_last_only){
    iter--; // To get the correct iteration number
    sprintf(buffer,"L");
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                        &str_end_idx);
  }

  free(L_phi);
  free(Cij);
  free(A);
  free(f);
  free(u);
  free(dx2);
  free(py1);
  free(py2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}

/*
- gigiaero, 07/04/2026, 1558 hours
*/
void solve_sor_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                     double *y,sim_parameters *config,
                                     bi_air_phys_mesh *b_a_m){
  // Solver variables
  double (*L_phi)[n-2] = calloc(m-2,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double *dx2 = malloc(sizeof(double)*(n-2));
  double *dy2 = malloc(sizeof(double)*(m-2));
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double *res = calloc(config->max_iter + 1,sizeof(double));
  FILE *file_log;

  // Configure log file
  sprintf(filename_log,"%s.log",config->casename);
  file_log = fopen(filename_log,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save,"%s_iter_",config->casename);
  find_str_end(filename_save,&str_end_idx);

  for(int i=1;i<n-1;i++){
    dx2[i-1] = delta_xy(x,i);
    dx2[i-1] *= dx2[i-1];
  }

  for(int j=1;j<m-1;j++){
    dy2[j-1] = delta_xy(y,j);
    dy2[j-1] *= dy2[j-1];
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res[iter])
          res[iter] = fabs(L_phi[j-1][i-1]);
      }
    }
    
    printf("SOR Iteration %010d | Res %.6e\n",iter,res[iter]);

    fprintf(file_log,"%.6e\n",res[iter]);

    // Test for convergence
    if(res[iter] <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res[iter] >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        Cij[j][i] = (-L_phi[j-1][i-1] - Cij[j][i-1]/dx2[i-1]
                     - Cij[j-1][i]/dy2[j-1])*config->r/
                     (-2./dx2[i-1] - 2./dy2[j-1]);
    }
    
    // Calculate phi^{n+1}
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++)
        phi[j][i] += Cij[j][i];
    }
    
    // Apply boundary conditions (low edge)
    apply_down_b_c(m,n,phi,x,y,b_a_m);

    if(!config->save_last_only)
      save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                          &str_end_idx);
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0 || config->save_last_only){
    iter--; // To get the correct iteration number
    sprintf(buffer,"L");
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,
                        &str_end_idx);
  }

  free(L_phi);
  free(Cij);
  free(dx2);
  free(dy2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  free(res);
  fclose(file_log);
}