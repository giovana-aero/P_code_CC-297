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
  double (*af_shape_dx)(double,double);

  if(b_a_m->af_type == 1)
    af_shape_dx = bi_air_shape_dx;

  if(b_a_m->af_type == 2)
    af_shape_dx = naca4_symm_dx;

  // if(b_a_m->af_type == 3)
  //   af_shape_dx = naca4_symm_dx_CST_t005;

  for(int i=1;i<n-1;i++){
    if(i >= b_a_m->ILE-1 && i < b_a_m->ITE){  // Airfoil
      if(i == b_a_m->ILE-1 && b_a_m->af_type == 2)
        phi[0][i] = phi[1][i] - (y[1] - y[0])*b_a_m->uinf*(*af_shape_dx)(b_a_m->t,x[i+1]*.5);
      else
        phi[0][i] = phi[1][i] - (y[1] - y[0])*b_a_m->uinf*(*af_shape_dx)(b_a_m->t,x[i]);
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

  // if(b_a_m->af_type == 2)
  //   x[b_a_m->ILE-1] += delta_x*.1;
}

/*
- gigiaero, 06/04/2026, 1051 hours
*/
double bi_air_shape(double t,double x_i){
  return 2.*t*x_i*(1. - x_i);
}

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

// /*
// - gigiaero, 09/04/2026, 1034 hours
// */
// void build_linear_sys_matrix_cols(int m,int n,double A[m-2][m-2],double *py1,
//                                   double *py2,double *dx2){
//   for(int i=1;i<n-1;i++){
//     A[0][0] = -1./dx2[i-1] - 1./py1[0] - 1./py2[0];
//     A[0][1] = 1./py1[0];

//     for(int j=1;j<m-3;j++){
//       A[j][j-1] = 1./py2[j];
//       A[j][j] = (-1./dx2[i-1] - 1./py1[j] - 1./py2[j]);
//       A[j][j+1] = 1./py1[j];
//     }

//     A[m-3][m-4] = 1./py2[m-3];
//     A[m-3][m-3] = -1./dx2[i-1] - 1./py1[m-3] - 1./py2[m-3];
//   }
// }

/*
this version calculates everything but the main diagonal. for lgs and slor only
- gigiaero, 09/04/2026, 1318 hours
*/
void build_linear_sys_matrix_cols2(int m,double A[m-2][m-2],double *py1,
                                  double *py2){
  // for(int i=1;i<n-1;i++){
  A[0][1] = 1./py1[0];

  for(int j=1;j<m-3;j++){
    A[j][j-1] = 1./py2[j];
    A[j][j+1] = 1./py1[j];
  }

  A[m-3][m-4] = 1./py2[m-3];
  // }
}

// /*
// for lgs and slor only
// - gigiaero, 10/04/2026, 1404 hours
// */
// void build_linear_sys_matrix_lines2(int m,int n,double A[n-2][n-2],double *px1,
//                                     double *px2){
//   for(int j=1;j<m-1;j++){
//     A[0][1] = 1./px1[0];

//     for(int i=1;i<n-3;i++){
//       A[i][i-1] = 1./px2[i];
//       A[i][i+1] = 1./px1[i];
//     }

//     A[n-3][n-4] = 1./px2[n-3];
//   }
// }

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

    case 6:
      puts("ADI, biconvex airfoil");
      // solve_adi_2d_rectangular_bi_air(m,n,phi,x,y,config,b_a_m);
      puts("implementation incomplete! :/ ");
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
(kudos to gibson for helping me find a certain dumb mistake)
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
- gigiaero, 10/04/2026, 2249 hours
*/
double naca4_symm_dx(double t,double x_i){
  // if(x_i == 0)
    // x_i = 1e-4;
    // return 5.*t*(-0.1260 - 2.*0.3516*x_i + 
          //  3.*0.2843*pow(x_i,2.) - 4.*0.1036*pow(x_i,3.));
    // return 0.;
  
  // else
    return 5.*t*(0.5*0.2969*pow(x_i,-0.5) - 0.1260 - 2.*0.3516*x_i + 
          3.*0.2843*pow(x_i,2.) - 4.*0.1036*pow(x_i,3.));
  // return 5.*t*(0.2969*pow(x_i,0.5) - 0.1260*x_i - 0.3516*pow(x_i,2.) + 
  //        0.2843*pow(x_i,3.) - 0.1036*pow(x_i,4.)); // original
}

/*
- gigiaero, 16/04/2026, 1940 hours
*/
double naca4_symm_CST_t005(double t,double x_i){
  // {0.002579,0.063728,0.070283,0.055543,0.064390,0.056486,
  //   3.523900,0}; original beta value in degrees
  double v_ex[] = {0.002579,0.063728,0.070283,0.055543,0.064390,0.056486,
                   0.061504,0};

  double A[] = {pow(2.*v_ex[0],.5),v_ex[1],v_ex[2],v_ex[3],v_ex[4],v_ex[5],
                tan(v_ex[6]),v_ex[7]};
  double N1 = .5,N2 = 1.;
  double sum = 0.,K;
  int n = 6;

  for(int r=0;r<=n;r++){
    K = factorial(n)/(factorial(r)*factorial(n-r));
    sum += A[r]*K*pow(x_i,r)*pow(1. - x_i,n - r);
  }

  return pow(x_i,N1)*pow(1.-x_i,N2)*sum;
}

/*
- gigiaero, 16/04/2026, 1940 hours
*/
double naca4_symm_dx_CST_t005(double t,double x_i){
  // {0.002579,0.063728,0.070283,0.055543,0.064390,0.056486,
  //   3.523900,0}; original beta value in degrees
  double v_ex[] = {0.002579,0.063728,0.070283,0.055543,0.064390,0.056486,
                   0.061504,0};

  double A[] = {pow(2.*v_ex[0],.5),v_ex[1],v_ex[2],v_ex[3],v_ex[4],v_ex[5],
                tan(v_ex[6]),v_ex[7]};
  double N1 = .5,N2 = 1.;
  double sum = 0.,K;
  int n = 6;

  for(int r=0;r<=n;r++){
    K = factorial(n)/(factorial(r)*factorial(n-r));
    sum += -A[r]*K*r*pow(x_i,r-1)*(n - r)*pow(1. - x_i,n - r - 1);
  }

  return -N1*pow(x_i,N1-1.)*N2*pow(1.-x_i,N2-1.)*sum;
}

/*
- gigiaero, 09/04/2026, 0857 hours
*/
void set_mesh_prmtrs(int mtype,bi_air_phys_mesh *b_a_mesh){
  switch(mtype){
    // ----------------------------------------------- original configuration
    case 1:
      b_a_mesh->ILE = 11;   // Leading edge
      b_a_mesh->ITE = 31;   // Trailing edge
      b_a_mesh->IMAX = 41;  // Number of points along x
      b_a_mesh->JMAX = 12;  // Number of points along y
      b_a_mesh->XSF = 1.25; // Stretching factor, x
      b_a_mesh->YSF = 1.25; // Stretching factor, y
      break;
    // ----------------------------------------------- 1/2
    case 2:
      b_a_mesh->ILE = 21;
      b_a_mesh->ITE = 61;
      b_a_mesh->IMAX = 81;
      b_a_mesh->JMAX = 24;
      b_a_mesh->XSF = 1.12237803;
      b_a_mesh->YSF =  1.10456837;
      break;
    // ----------------------------------------------- 1/4
    case 3: 
      b_a_mesh->ILE = 41;
      b_a_mesh->ITE = 121;
      b_a_mesh->IMAX = 161;
      b_a_mesh->JMAX = 48;
      b_a_mesh->XSF = 1.0605144;
      b_a_mesh->YSF = 1.04818314;
      break;
    // ----------------------------------------------- 1/8
    case 4: 
      b_a_mesh->ILE = 81;
      b_a_mesh->ITE = 241;
      b_a_mesh->IMAX = 321;
      b_a_mesh->JMAX = 96;
      b_a_mesh->XSF = 1.03008612;
      b_a_mesh->YSF = 1.0231651;
      break;
    // ----------------------------------------------- 1/.5
    case 5:
      b_a_mesh->ILE = 6;
      b_a_mesh->ITE = 16;
      b_a_mesh->IMAX = 21;
      b_a_mesh->JMAX = 6;
      b_a_mesh->XSF = 1.51979153;
      b_a_mesh->YSF = 1.77743291;
      break;

    default:
      puts("Invalid mtype");
      exit(4);
  }
}

void solve_adi_2d_rectangular_bi_air(int m,int n,double phi[m][n],double *x,
                                     double *y,sim_parameters *config,
                                     bi_air_phys_mesh *b_a_m){
  // Solver variables
  double (*L_phi)[n-2] = calloc(m-2,sizeof *L_phi);
  double (*Cij)[n] = calloc(m,sizeof *Cij);
  double (*fij)[n] = calloc(m,sizeof *fij);
  double (*A1)[n-2] = calloc(n-2,sizeof *A1);
  double (*A2)[m-2] = calloc(m-2,sizeof *A2);
  double *f1 = malloc(sizeof(double)*(n-2));
  double *u1 = malloc(sizeof(double)*(n-2));
  double *f2 = malloc(sizeof(double)*(m-2));
  double *u2 = malloc(sizeof(double)*(m-2));
  double *px1 = malloc(sizeof(double)*(n-2));
  double *px2 = malloc(sizeof(double)*(n-2));
  double *py1 = malloc(sizeof(double)*(m-2));
  double *py2 = malloc(sizeof(double)*(m-2));
  double vars_old,val;
  double a = config->r;
  double w = config->w;
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double res;
  FILE *file_log;

  // Configure log file
  sprintf(filename_log,"%s.log",config->casename);
  file_log = fopen(filename_log,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save,"%s_iter_",config->casename);
  find_str_end(filename_save,&str_end_idx);

  // for(int i=1;i<n-1;i++){
  //   dx2[i-1] = delta_xy(x,i);
  //   dx2[i-1] *= dx2[i-1];
  // }

  for(int i=1;i<n-1;i++){
    px1[i-1] = (x[i+1] - x[i-1])*(x[i+1] - x[i]);
    px2[i-1] = (x[i+1] - x[i-1])*(x[i] - x[i-1]);
  }

  for(int j=1;j<m-1;j++){
    py1[j-1] = (y[j+1] - y[j-1])*(y[j+1] - y[j]);
    py2[j-1] = (y[j+1] - y[j-1])*(y[j] - y[j-1]);
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  // for(int j=1;j<n-1;j++){
  //   for(int i=1;i<n-1;i++)
  //       fij[j][i] = 1.;
  // }

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    res = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res)
          res = fabs(L_phi[j-1][i-1]);
      }
    }

    // print_2d_array(m-2,n-2,L_phi);

    printf("ADI Iteration %010d | Res %.6e\n",iter,res);

    fprintf(file_log,"%.6e\n",res);

    // Test for convergence
    if(res <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for fij
    for(int j=1;j<m-1;j++){
      A1[0][0] =  a + 2./px1[0] + 2./px2[0];
      A1[0][1] = -2./px1[0];
      f1[0] = a*w*L_phi[j][0] + 2./px2[0]*fij[j][0];
      // printf("%f\n",L_phi[j][1]);
  
      for(int i=1;i<n-3;i++){
        A1[i][i-1] = -2./px2[i];
        A1[i][i] = a + 2./px1[i] + 2./px2[i];
        A1[i][i+1] = -2./px1[i];
        f1[i] = a*w*L_phi[j][i];
      }
  
      A1[n-3][n-4] = -2./px2[n-3];
      A1[n-3][n-3] = a + 2./px1[n-3] + 2./px2[n-3];
      f1[n-3] = a*w*L_phi[j][n-3] + 2./px1[n-3]*fij[j][n-1];

      diagonal_matrix_solver(n-2,A1,f1,u1);

      // print_1d_array(n-2,u1);

      for(int i=1;i<n-1;i++)
        fij[j][i] = u1[i-1];

    }

    // // print_2d_array(n-2,n-2,A1);
    // print_1d_array(n-2,f1);
    // print_1d_array(n-2,u1);
    // putchar('\n');

    // print_2d_array_to_file(n-2,n-2,A1,"A1.txt",0);

    // Solve for Cij
    for(int i=1;i<n-1;i++){
      A2[0][0] = a + 2./py1[0] + 2./py2[0];
      A2[0][1] = -2./py1[0];
      f2[0] = fij[1][i] + 2./py2[0]*Cij[1][i]; // implicit bc on Cij

      for(int j=1;j<m-3;j++){
        A2[j][j-1] = -2./py2[j];
        A2[j][j] = a + 2./py1[j] + 2./py2[j];
        A2[j][j+1] = -2./py1[j];
        f2[j] = fij[j+1][i];
      }
      
      A2[m-3][m-4] = -2./py2[m-3];
      A2[m-3][m-3] = a + 2./py1[m-3] + 2./py2[m-3];
      f2[m-3] = fij[m-3][i] + 2./py1[m-3]*Cij[m-1][i];
      
      diagonal_matrix_solver(m-2,A2,f2,u2);

      for(int j=1;j<m-1;j++)
        Cij[j][i] = u2[j-1];
    }

    // // print_2d_array(m-2,m-2,A2);
    // print_1d_array(m-2,f2);
    // print_1d_array(m-2,u2);
    // // print_2d_array_to_file(m-2,m-2,A2,"A2.txt",0);

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
  free(fij);
  free(A1);
  free(A2);
  free(f1);
  free(u1);
  free(f2);
  free(u2);
  free(px1);
  free(px2);
  free(py1);
  free(py2);
  free(filename_save);
  free(filename_log);
  free(buffer);
  fclose(file_log);
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
  double res;
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
    res = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res)
          res = fabs(L_phi[j-1][i-1]);
      }
    }
    
    printf("GS Iteration %010d | Res %.6e\n",iter,res);

    fprintf(file_log,"%.6e\n",res);

    // Test for convergence
    if(res <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int i=1;i<n-1;i++){
      Cij[1][i] = (-L_phi[0][i-1] - Cij[1][i-1]/dx2[i-1] - Cij[1][i]/dy2[0])/
                  (-2./dx2[i-1] - 2./dy2[0]);
      for(int j=2;j<m-1;j++)
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
  double res;
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
    py1[j-1] = (y[j+1] - y[j-1])*(y[j+1] - y[j]);
    py2[j-1] = (y[j+1] - y[j-1])*(y[j] - y[j-1]);
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    res = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res)
          res = fabs(L_phi[j-1][i-1]);
      }
    }

    printf("LGS Iteration %010d | Res %.6e\n",iter,res);

    fprintf(file_log,"%.6e\n",res);

    // Test for convergence
    if(res <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Build the matrix of the linear system
    build_linear_sys_matrix_cols2(m,A,py1,py2);

    // Solve for Cij
    for(int i=1;i<n-1;i++){
      A[0][0] = -1./dx2[i-1] - 1./py1[0] - 1./py2[0];               
      f[0] = (-L_phi[0][i-1] - Cij[1][i]/dx2[i-1])/2. - Cij[1][i]/py2[0];

      for(int j=1;j<m-3;j++){
        A[j][j] = (-1./dx2[i-1] - 1./py1[j] - 1./py2[j]);
        f[j] = (-L_phi[j][i-1] - Cij[j+1][i]/dx2[i-1])/2.;
      }

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
  int iter = 0;
  // Save files
  char *filename_save = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log = malloc(sizeof(char)*200);
  double res;
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
    res = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j][i],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j][i]) > res)
          res = L_phi[j][i];
      }
    }
    
    printf("PJ Iteration %010d | Res %.6e\n",iter,res);

    fprintf(file_log,"%.6e\n",res);

    // Test for convergence
    if(res <= config->eps & iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res >= div_ref){
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
  double res;
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
    py1[j-1] = (y[j+1] - y[j-1])*(y[j+1] - y[j]);
    py2[j-1] = (y[j+1] - y[j-1])*(y[j] - y[j-1]);
  }

  // Save initial condition
  if(config->save_i_c)
    save_results_qtimes(m,n,phi,&iter,config,buffer,filename_save,&str_end_idx);

  build_linear_sys_matrix_cols2(m,A,py1,py2);

  for(iter;iter<=config->max_iter;iter++){
    // Calculate residual operator
    res = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res)
          res = fabs(L_phi[j-1][i-1]);
      }
    }

    printf("SLOR Iteration %010d | Res %.6e\n",iter,res);

    fprintf(file_log,"%.6e\n",res);

    // Test for convergence
    if(res <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int i=1;i<n-1;i++){
      A[0][0] = -1./dx2[i-1] - 1./py1[0] - 1./py2[0];               
      f[0] = (-L_phi[0][i-1] - Cij[1][i]/dx2[i-1])*r/2. - Cij[1][i]/py2[0];

      for(int j=1;j<m-3;j++){
        A[j][j] = (-1./dx2[i-1] - 1./py1[j] - 1./py2[j]);
        f[j] = (-L_phi[j][i-1] - Cij[j+1][i]/dx2[i-1])*r/2.;
      }

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
  double res;
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
    res = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        scheme_der2_o2_central_var_deltas_xy(&L_phi[j-1][i-1],m,n,phi,x,y,i,j);
        
        if(fabs(L_phi[j-1][i-1]) > res)
          res = fabs(L_phi[j-1][i-1]);
      }
    }
    
    printf("SOR Iteration %010d | Res %.6e\n",iter,res);

    fprintf(file_log,"%.6e\n",res);

    // Test for convergence
    if(res <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for Cij
    for(int i=1;i<n-1;i++){
      Cij[1][i] = (-L_phi[0][i-1] - Cij[1][i-1]/dx2[i-1] - Cij[1][i]/dy2[0])*
                   config->r/(-2./dx2[i-1] - 2./dy2[0]);
      for(int j=2;j<m-1;j++)
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
  fclose(file_log);
}