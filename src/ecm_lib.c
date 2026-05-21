#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/num_methods.h"
#include"../include/ecm_lib.h"
#include"../include/eom_lib.h"
#include"../include/strings.h"

#define div_ref 1e100
#define pi 3.1415926535897932384626433


void calc_A_ecm(int m,int n,double A[m][n],double x[m][n],double y[m][n],
                int TE){
  for(int i=1;i<TE-1;i++)
    A[0][i] = pow(uniform_scheme_der1_o2_central_ecm(m,n,x,i),2.) + 
                pow(uniform_scheme_der1_o2_central_ecm(m,n,y,i),2.);

  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++)
      A[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,2),2.) + 
                pow(uniform_scheme_der1_o2_central(m,n,y,i,j,2),2.);
  }
}

/*
- gigiaero, 20/05/2026, 1457 hours
*/
void calc_B_ecm(int m,int n,double B[m][n],double x[m][n],double y[m][n],
                int TE){
  for(int i=1;i<TE-1;i++)                
  B[0][i] = uniform_scheme_der1_o2_central(m,n,x,i,0,1)*
            uniform_scheme_der1_o2_central_ecm(m,n,x,i) + 
            uniform_scheme_der1_o2_central(m,n,y,i,0,1)*
            uniform_scheme_der1_o2_central_ecm(m,n,y,i);
  
  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++)
      B[j][i] = uniform_scheme_der1_o2_central(m,n,x,i,j,1)*
                uniform_scheme_der1_o2_central(m,n,x,i,j,2) + 
                uniform_scheme_der1_o2_central(m,n,y,i,j,1)*
                uniform_scheme_der1_o2_central(m,n,y,i,j,2);
  }
}

/*
- gigiaero, 20/05/2026, 1441 hours
*/
void calc_C_ecm(int m,int n,double C[m][n],double x[m][n],double y[m][n],
                int TE){
  for(int i=1;i<TE-1;i++)
    C[0][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,0,1),2.) + 
              pow(uniform_scheme_der1_o2_central(m,n,y,i,0,1),2.);

  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++)
      C[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,1),2.) + 
                pow(uniform_scheme_der1_o2_central(m,n,y,i,j,1),2.);
  }
}

/*
- gigiaero, 20/05/2026, 1438 hours
*/
void evaluate_delta_form_ecm(sim_prmtrs *config,msh_prmtrs *msh,
                             control_prmtrs *c_prmtrs,int init_only){
  int n = msh->IMAX+2*msh->JMAX;
  int m = msh->JMAX;
  double (*x)[n] = calloc(m,sizeof *x);
  double (*y)[n] = calloc(m,sizeof *y);

  initialize_mesh_ecm(m,n,x,y,msh);

  if(config->save_i_c){
    char *filename = malloc(sizeof(char)*200);

    sprintf(filename,"%s%s",config->casename,"_x_initial.dat");
    print_2d_array_to_file(m,n,x,filename,0);
    sprintf(filename,"%s%s",config->casename,"_y_initial.dat");
    print_2d_array_to_file(m,n,y,filename,0);

    free(filename);
  }

  if(!init_only){
    puts("ADI, elliptical C mesh");
    solve_adi_2d_rectangular_ecm(m,n,x,y,config,msh->JMAX+1);
  }

  free(x);
  free(y);
}

/*
- gigiaero, 20/05/2026, 1438 hours
*/
void half_ellipse(double *x,double *y,double *prmtrs,int n,int invert_th){
  double *th = malloc(sizeof(double)*n);

  if(invert_th)
    linspace(th,3.*pi/2.,pi/2.,n);
  else
    linspace(th,pi/2.,3.*pi/2.,n);

  for(int i=0;i<n;i++){
    x[i] = prmtrs[2] + prmtrs[0]*cos(th[i]);
    y[i] = prmtrs[3] + prmtrs[1]*sin(th[i]);
  }
    
  free(th);
}

/*
- gigiaero, 20/05/2026, 1438 hours
*/
void initialize_mesh_ecm(int m,int n,double x[m][n],double y[m][n],
                       msh_prmtrs *msh){
  int chord_n = (msh->IMAX+1)/2;
  double *x_axis = malloc(sizeof(double)*chord_n);
  double *tmp_x = malloc(sizeof(double)*msh->IMAX);
  double *tmp_y = malloc(sizeof(double)*msh->IMAX);
  double *tmp_x_wake = malloc(sizeof(double)*(m+1));
  cosspace(x_axis,0.,msh->c,chord_n,0);

  msh->end_prmtrs[2] = msh->c;
  
  switch(msh->af_type){
    case 1:
      init_af_bi_air(tmp_x,tmp_y,x_axis,chord_n,msh);
      break;

    case 2:
      init_af_naca4(tmp_x,tmp_y,x_axis,chord_n,msh);
      break;

    case 3:
      init_af_cst(tmp_x,tmp_y,x_axis,chord_n,msh);
      break;

    default:
      puts("initialize_mesh: invalid af_type");
      exit(13);
  }

  for(int i=0,k=m;i<msh->IMAX;i++,k++){
    x[0][k] = tmp_x[i];
    y[0][k] = tmp_y[i];
    x[m-1][k] = tmp_x[i];
    y[m-1][k] = tmp_y[i];
  }

  half_ellipse(tmp_x,tmp_y,msh->end_prmtrs,msh->IMAX,1);

  for(int i=0,k=m;i<msh->IMAX;i++,k++){
    x[m-1][k] = tmp_x[i];
    y[m-1][k] = tmp_y[i];
  }

  cosspace(tmp_x_wake,msh->end_prmtrs[2],msh->c+msh->end_prmtrs[0],m+1,1);
  
  for(int i=0,k1=m,k2=msh->IMAX+m;i<m;i++,k1--,k2++){
    x[0][i] = tmp_x_wake[k1];
    x[0][k2] = tmp_x_wake[i+1];
    y[0][i] = 0.;
    y[0][k2] = 0.;
    x[m-1][i] = tmp_x_wake[k1];
    x[m-1][k2] = tmp_x_wake[i+1];
    y[m-1][i] = -msh->end_prmtrs[0];
    y[m-1][k2] = msh->end_prmtrs[0];
  }

  init_type3(m,n,x,y,msh);

  free(x_axis);
  free(tmp_x);
  free(tmp_y);
  free(tmp_x_wake);
}

/*
- gigiaero, 20/05/2026, 1502 hours
*/
void L_phi_ecm(int m,int n,double L_phi_x[m][n],double L_phi_y[m][n],
               double x[m][n],double y[m][n],double A[m][n],
               double B[m][n],double C[m][n],int TE){
  for(int i=1;i<TE-1;i++){
    L_phi_x[0][i] = A[0][i]*uniform_scheme_der2_o2_central(m,n,x,i,0,1) - 
                    2.*B[0][i]*uniform_scheme_der2_o2_central_ecm(m,n,x,i,3) + 
                    C[0][i]*uniform_scheme_der2_o2_central_ecm(m,n,x,i,2);
    // disp(A[0][i]);
    // disp(B[0][i]);
    // disp(C[0][i]);
    // disp(uniform_scheme_der2_o2_central(m,n,x,i,0,1));
    // disp(uniform_scheme_der2_o2_central_ecm(m,n,x,i,3));
    // disp(uniform_scheme_der2_o2_central_ecm(m,n,x,i,2));
    L_phi_y[0][i] = A[0][i]*uniform_scheme_der2_o2_central(m,n,y,i,0,1) - 
                    2.*B[0][i]*uniform_scheme_der2_o2_central_ecm(m,n,y,i,3) + 
                    C[0][i]*uniform_scheme_der2_o2_central_ecm(m,n,y,i,2);  
  }

  for(int j=1;j<m-1;j++){
    for(int i=1;i<n-1;i++){
      L_phi_x[j][i] = A[j][i]*uniform_scheme_der2_o2_central(m,n,x,i,j,1) - 
                       2.*B[j][i]*uniform_scheme_der2_o2_central(m,n,x,i,j,3) + 
                       C[j][i]*uniform_scheme_der2_o2_central(m,n,x,i,j,2);
                      //  D[j][i]*(P*uniform_scheme_der1_o2_central(m,n,x,i,j,1) +
                      //           Q*uniform_scheme_der1_o2_central(m,n,x,i,j,2));
                       
      L_phi_y[j][i] = A[j][i]*uniform_scheme_der2_o2_central(m,n,y,i,j,1) - 
                       2.*B[j][i]*uniform_scheme_der2_o2_central(m,n,y,i,j,3) + 
                       C[j][i]*uniform_scheme_der2_o2_central(m,n,y,i,j,2);
                      //  D[j][i]*(P*uniform_scheme_der1_o2_central(m,n,y,i,j,1) +
                      //           Q*uniform_scheme_der1_o2_central(m,n,y,i,j,2));
    }
  }
}

void solve_adi_2d_rectangular_ecm(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config,int TE){
  // Solver variables
  int m_wake = (2*m-1);
  double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
  double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
  double (*Delta_x)[n] = calloc(m,sizeof *Delta_x);
  double (*Delta_y)[n] = calloc(m,sizeof *Delta_y);
  double (*A)[n] = calloc(m,sizeof *A);
  double (*B)[n] = calloc(m,sizeof *B);
  double (*C)[n] = calloc(m,sizeof *C);
  // double (*D)[n] = calloc(m,sizeof *D);
  double (*fx)[n] = calloc(m,sizeof *fx);
  double (*fy)[n] = calloc(m,sizeof *fy);
  double *ak = malloc(sizeof(double)*(n-3));
  double *bk = malloc(sizeof(double)*(n-2));
  double *ck = malloc(sizeof(double)*(n-3));
  double *fkx = malloc(sizeof(double)*(n-2));
  double *fky = malloc(sizeof(double)*(n-2));
  double *ukx = malloc(sizeof(double)*(n-2));
  double *uky = malloc(sizeof(double)*(n-2));
  double *an = malloc(sizeof(double)*(m_wake-3));
  double *bn = malloc(sizeof(double)*(m_wake-2));
  double *cn = malloc(sizeof(double)*(m_wake)-3);
  double *fnx = malloc(sizeof(double)*(m_wake-2));
  double *fny = malloc(sizeof(double)*(m_wake-2));
  double *unx = malloc(sizeof(double)*(m_wake-2));
  double *uny = malloc(sizeof(double)*(m_wake-2));
  double omega = config->w;
  double alpha = config->alpha;
  int k = 1;
  int iter = 1;
  // Save files
  char *filename_save_x = malloc(sizeof(char)*200);
  char *filename_save_y = malloc(sizeof(char)*200);
  char *buffer = malloc(sizeof(char)*200);
  int str_end_idx;
  // Residuals
  char *filename_log_x = malloc(sizeof(char)*200);
  char *filename_log_y = malloc(sizeof(char)*200);
  double res_x,res_y;
  FILE *file_log_x;
  FILE *file_log_y;
   
  // Configure log files
  sprintf(filename_log_x,"%s_x.log",config->casename);
  file_log_x = fopen(filename_log_x,"w");
  sprintf(filename_log_y,"%s_y.log",config->casename);
  file_log_y = fopen(filename_log_y,"w");

  // Prepare string to save simulation data  
  sprintf(filename_save_x,"%s_x_iter_",config->casename);
  sprintf(filename_save_y,"%s_y_iter_",config->casename);
  find_str_end(filename_save_x,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    calc_A_ecm(m,n,A,x,y,TE);
    calc_B_ecm(m,n,B,x,y,TE);
    calc_C_ecm(m,n,C,x,y,TE);

    // print_2d_array_to_file(m,n,A,"./mat_A.dat",1);
    // print_2d_array_to_file(m,n,B,"./mat_B.dat",1);
    // print_2d_array_to_file(m,n,C,"./mat_C.dat",1);

    L_phi_ecm(m,n,L_phi_x,L_phi_y,x,y,A,B,C,TE);

    // print_2d_array_to_file(m,n,L_phi_x,"./mat_L_phi_x.dat",1);
    // print_2d_array_to_file(m,n,L_phi_y,"./mat_L_phi_y.dat",1);

    res_x = 0.;
    res_y = 0.;
    for(int i=1;i<TE;i++){
      if(fabs(L_phi_x[0][i]) > res_x)
        res_x = fabs(L_phi_x[0][i]);

      if(fabs(L_phi_y[0][i]) > res_y)
        res_y = fabs(L_phi_y[0][i]);
    }
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        if(fabs(L_phi_x[j][i]) > res_x)
          res_x = fabs(L_phi_x[j][i]);

        if(fabs(L_phi_y[j][i]) > res_y)
          res_y = fabs(L_phi_y[j][i]);
      }
    }

    printf("ADI Iteration %010d | Res x %.6e | Res y %.6e\n",iter,res_x,res_y);

    fprintf(file_log_x,"%.6e\n",res_x);
    fprintf(file_log_y,"%.6e\n",res_y);

    // Test for convergence
    if(res_x <= config->eps && res_y <= config->eps && iter != 0){
      puts("<< Convergence! >>");
      iter++;
      break;
    }

    if(res_x >= div_ref || res_y >= div_ref){
      puts("- Divergence");
      iter++;
      break;
    }

    // Solve for the deltas - step 1 (eta)
    for(int i1=1,i2=n-2,k;i1<TE;i1++,i2--){
      // disp(n);
      // disp(i1);
      // disp(i2);
      k = 0;
      // printf("%d | %f %f\n",i,x[0][i],y[0][i]);
      
      bn[k] = alpha + 2.*C[m-2][i2];
      cn[k] = -C[m-2][i2];

      fnx[k] = alpha*omega*L_phi_x[m-2][i2];
      fny[k] = alpha*omega*L_phi_y[m-2][i2];
      
      k++;

      for(int j=m-3;j>=1;j--,k++){
        // disp(k);
        an[k-1] = -C[j][i2];
        bn[k] = alpha + 2.*C[j][i2];
        cn[k] = -C[j][i2];

        fnx[k] = alpha*omega*L_phi_x[j][i2];
        fny[k] = alpha*omega*L_phi_y[j][i2];
      }

      for(int j=0;j<m-2;j++,k++){
        // disp(k);
        an[k-1] = -C[j][i1];
        bn[k] = alpha + 2.*C[j][i1];
        cn[k] = -C[j][i1];

        fnx[k] = alpha*omega*L_phi_x[j][i1];
        fny[k] = alpha*omega*L_phi_y[j][i1];
      }
      // disp(m_wake-3);

      an[m_wake-4] = -C[m-2][i1];
      bn[m_wake-3] = alpha + 2.*C[m-2][i1];

      fnx[m_wake-3] = alpha*omega*L_phi_x[m-2][i1];
      fny[m_wake-3] = alpha*omega*L_phi_y[m-2][i1];

      tridiagonal_matrix_solver(m_wake-2,an,bn,cn,fnx,unx);
      tridiagonal_matrix_solver(m_wake-2,an,bn,cn,fny,uny);

      for(int j=1,k1=(m_wake-1)/2,k2=k1-2;j<m-1;j++,k1++,k2--){
        // printf("%02d %02d %02d\n",j,k1,k2);
        fx[j][i1] = unx[k1];
        fy[j][i1] = uny[k1];
        fx[j][i2] = unx[k2];
        fy[j][i2] = uny[k2];
      }

      // disp((m_wake-3)/2);
      // disp((m_wake-1)/2);
      // return;
      fx[0][i1] = unx[(m_wake-3)/2];
      fy[0][i1] = uny[(m_wake-3)/2];
      fx[0][i2] = unx[(m_wake-3)/2];
      fy[0][i2] = uny[(m_wake-3)/2];
    }

    for(int i=TE;i<n-TE;i++){
      bn[0] = alpha + 2.*C[1][i];
      cn[0] = -C[1][i];

      fnx[0] = alpha*omega*L_phi_x[1][i];
      fny[0] = alpha*omega*L_phi_y[1][i];
      
      for(int j=2;j<m-2;j++){
        an[j-2] = -C[j][i];
        bn[j-1] = alpha + 2.*C[j][i];
        cn[j-1] = -C[j][i];

        fnx[j-1] = alpha*omega*L_phi_x[j][i];
        fny[j-1] = alpha*omega*L_phi_y[j][i];
      }

      an[m-4] = -C[m-2][i];
      bn[m-3] = alpha + 2.*C[m-2][i];

      fnx[m-3] = alpha*omega*L_phi_x[m-2][i];
      fny[m-3] = alpha*omega*L_phi_y[m-2][i];

      tridiagonal_matrix_solver(m-2,an,bn,cn,fnx,unx);
      tridiagonal_matrix_solver(m-2,an,bn,cn,fny,uny);

      for(int j=1;j<m-1;j++){
        fx[j][i] = unx[j-1];
        fy[j][i] = uny[j-1];
      }
    }

    // // Reapply "periodicity" boundary condition
    // for(int i1=1,i2=n-2;i1<TE;i1++,i2--){
    //   fx[0][i2] = fx[0][i1];
    //   fy[0][i2] = fy[0][i1];
    // }

    // Solve for the deltas - step 2 (ksi)
    bk[0] = alpha + 2.*A[0][1];
    ck[0] = -A[0][1];

    fkx[0] = fx[0][1];
    fky[0] = fy[0][1];
    
    for(int i=2;i<TE-2;i++){
      ak[i-2] = -A[0][i];
      bk[i-1] = alpha + 2.*A[0][i];
      ck[i-1] = -A[0][i];

      fkx[i-1] = fx[0][i];
      fky[i-1] = fy[0][i];
    }

    ak[TE-4] = -A[0][TE-2];
    bk[TE-3] = alpha + 2.*A[0][TE-2];

    fkx[TE-3] = fx[0][TE-2];
    fky[TE-3] = fy[0][TE-2];
    
    tridiagonal_matrix_solver(TE-2,ak,bk,ck,fkx,ukx);
    tridiagonal_matrix_solver(TE-2,ak,bk,ck,fky,uky);

    // print_1d_array(TE-4,ak);
    // print_1d_array(TE-3,bk);
    // print_1d_array(TE-4,ck);
    // print_1d_array(TE-3,fkx);
    // print_1d_array(TE-3,fky);
    // print_1d_array(TE-3,ukx);
    // print_1d_array(TE-3,uky);

    // disp(TE);return;
    for(int i=1;i<TE-1;i++){
      Delta_x[0][i] = ukx[i-1];
      Delta_y[0][i] = uky[i-1];
    }

    // for(int i=1;i<TE-1;i++){
    //   Delta_x[0][i] = fx[0][i]/alpha;
    //   Delta_y[0][i] = fy[0][i]/alpha;
    // }

    for(int j=1;j<m-1;j++){
      bk[0] = alpha + 2.*A[j][1];
      ck[0] = -A[j][1];

      fkx[0] = fx[j][1];
      fky[0] = fy[j][1];
      
      for(int i=2;i<n-2;i++){
        ak[i-2] = -A[j][i];
        bk[i-1] = alpha + 2.*A[j][i];
        ck[i-1] = -A[j][i];

        fkx[i-1] = fx[j][i];
        fky[i-1] = fy[j][i];
      }

      ak[n-4] = -A[j][n-2];
      bk[n-3] = alpha + 2.*A[j][n-2];

      fkx[n-3] = fx[j][n-2];
      fky[n-3] = fy[j][n-2];
      
      tridiagonal_matrix_solver(n-2,ak,bk,ck,fkx,ukx);
      tridiagonal_matrix_solver(n-2,ak,bk,ck,fky,uky);

      for(int i=1;i<n-1;i++){
        Delta_x[j][i] = ukx[i-1];
        Delta_y[j][i] = uky[i-1];
      }
    }

    // Calculate new x and y
    for(int i=1;i<TE;i++){
      x[0][i] += Delta_x[0][i];
      y[0][i] += Delta_y[0][i];
    }

    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        x[j][i] += Delta_x[j][i];
        y[j][i] += Delta_y[j][i];
      }
    }

    // print_2d_array_to_file(m,n,A,"./mat_A.dat",1);
    // print_2d_array_to_file(m,n,B,"./mat_B.dat",1);
    // print_2d_array_to_file(m,n,C,"./mat_C.dat",1);
    // print_2d_array_to_file(m,n,L_phi_x,"./mat_L_phi_x.dat",1);
    // print_2d_array_to_file(m,n,L_phi_y,"./mat_L_phi_y.dat",1);
    // print_2d_array_to_file(m,n,fx,"./mat_fx.dat",1);
    // print_2d_array_to_file(m,n,fy,"./mat_fy.dat",1);
    // print_2d_array_to_file(m,n,Delta_x,"./mat_deltax.dat",1);
    // print_2d_array_to_file(m,n,Delta_y,"./mat_deltay.dat",1);
    // exit(1);

    // Reapply "periodicity" boundary condition
    for(int i1=1,i2=n-2;i1<TE;i1++,i2--){
      x[0][i2] = x[0][i1];
      y[0][i2] = y[0][i1];
    }

    if(!config->save_last_only){
      save_results_qtimes(m,n,x,&iter,config,buffer,filename_save_x,
                          &str_end_idx);
      save_results_qtimes(m,n,y,&iter,config,buffer,filename_save_y,
                          &str_end_idx);
    }
  }

  // Save last iteration if it wasn't saved
  if(iter%config->qtimes != 0 || config->save_last_only){
    iter--; // To get the correct iteration number
    sprintf(buffer,"L");
    save_results_qtimes(m,n,x,&iter,config,buffer,filename_save_x,&str_end_idx);
    sprintf(buffer,"L");
    save_results_qtimes(m,n,y,&iter,config,buffer,filename_save_y,&str_end_idx);
  }

  free(L_phi_x);
  free(L_phi_y);
  free(Delta_x);
  free(Delta_y);
  free(A);
  free(B);
  free(C);
  // free(D);
  free(fx);
  free(fy);
  free(ak);
  free(bk);
  free(ck);
  free(ukx);
  free(uky);
  free(fkx);
  free(fky);
  free(an);
  free(bn);
  free(cn);
  free(unx);
  free(uny);
  free(fnx);
  free(fny);
  free(filename_save_x);
  free(filename_save_y);
  free(buffer);
  free(filename_log_x);
  free(filename_log_y);
  fclose(file_log_x);
  fclose(file_log_y);
}

// /*
// step 1: ksi, step 2: eta
// (not working)
// - gigiaero, 21/05/2026, 1056 hours
// */
// void solve_adi_2d_rectangular_ecm(int m,int n,double x[m][n],double y[m][n],
//                                   sim_prmtrs *config,int TE){
//   // Solver variables
//   int m_wake = (2*m-1);
//   double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
//   double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
//   double (*Delta_x)[n] = calloc(m,sizeof *Delta_x);
//   double (*Delta_y)[n] = calloc(m,sizeof *Delta_y);
//   double (*A)[n] = calloc(m,sizeof *A);
//   double (*B)[n] = calloc(m,sizeof *B);
//   double (*C)[n] = calloc(m,sizeof *C);
//   // double (*D)[n] = calloc(m,sizeof *D);
//   double (*fx)[n] = calloc(m,sizeof *fx);
//   double (*fy)[n] = calloc(m,sizeof *fy);
//   double *ak = malloc(sizeof(double)*(n-3));
//   double *bk = malloc(sizeof(double)*(n-2));
//   double *ck = malloc(sizeof(double)*(n-3));
//   double *fkx = malloc(sizeof(double)*(n-2));
//   double *fky = malloc(sizeof(double)*(n-2));
//   double *ukx = malloc(sizeof(double)*(n-2));
//   double *uky = malloc(sizeof(double)*(n-2));
//   double *an = malloc(sizeof(double)*(m_wake-3));
//   double *bn = malloc(sizeof(double)*(m_wake-2));
//   double *cn = malloc(sizeof(double)*(m_wake)-3);
//   double *fnx = malloc(sizeof(double)*(m_wake-2));
//   double *fny = malloc(sizeof(double)*(m_wake-2));
//   double *unx = malloc(sizeof(double)*(m_wake-2));
//   double *uny = malloc(sizeof(double)*(m_wake-2));
//   double omega = config->w;
//   double alpha = config->alpha;
//   int k = 1;
//   int iter = 1;
//   // Save files
//   char *filename_save_x = malloc(sizeof(char)*200);
//   char *filename_save_y = malloc(sizeof(char)*200);
//   char *buffer = malloc(sizeof(char)*200);
//   int str_end_idx;
//   // Residuals
//   char *filename_log_x = malloc(sizeof(char)*200);
//   char *filename_log_y = malloc(sizeof(char)*200);
//   double res_x,res_y;
//   FILE *file_log_x;
//   FILE *file_log_y;
   
//   // Configure log files
//   sprintf(filename_log_x,"%s_x.log",config->casename);
//   file_log_x = fopen(filename_log_x,"w");
//   sprintf(filename_log_y,"%s_y.log",config->casename);
//   file_log_y = fopen(filename_log_y,"w");

//   // Prepare string to save simulation data  
//   sprintf(filename_save_x,"%s_x_iter_",config->casename);
//   sprintf(filename_save_y,"%s_y_iter_",config->casename);
//   find_str_end(filename_save_x,&str_end_idx);

//   for(iter;iter<=config->max_iter;iter++){
//     calc_A_ecm(m,n,A,x,y,TE);
//     calc_B_ecm(m,n,B,x,y,TE);
//     calc_C_ecm(m,n,C,x,y,TE);

//     print_2d_array_to_file(m,n,A,"./mat_A.dat",1);
//     print_2d_array_to_file(m,n,B,"./mat_B.dat",1);
//     print_2d_array_to_file(m,n,C,"./mat_C.dat",1);

//     L_phi_ecm(m,n,L_phi_x,L_phi_y,x,y,A,B,C,TE);

//     print_2d_array_to_file(m,n,L_phi_x,"./mat_L_phi_x.dat",1);
//     print_2d_array_to_file(m,n,L_phi_y,"./mat_L_phi_y.dat",1);

//     res_x = 0.;
//     res_y = 0.;
//     for(int i=0;i<TE;i++){
//       if(fabs(L_phi_x[0][i]) > res_x)
//         res_x = fabs(L_phi_x[0][i]);

//       if(fabs(L_phi_y[0][i]) > res_y)
//         res_y = fabs(L_phi_y[0][i]);
//     }
//     for(int j=1;j<m-1;j++){
//       for(int i=0;i<n-1;i++){
//         if(fabs(L_phi_x[j][i]) > res_x)
//           res_x = fabs(L_phi_x[j][i]);

//         if(fabs(L_phi_y[j][i]) > res_y)
//           res_y = fabs(L_phi_y[j][i]);
//       }
//     }

//     printf("ADI Iteration %010d | Res x %.6e | Res y %.6e\n",iter,res_x,res_y);

//     fprintf(file_log_x,"%.6e\n",res_x);
//     fprintf(file_log_y,"%.6e\n",res_y);

//     // Test for convergence
//     if(res_x <= config->eps && res_y <= config->eps && iter != 0){
//       puts("<< Convergence! >>");
//       iter++;
//       break;
//     }

//     if(res_x >= div_ref || res_y >= div_ref){
//       puts("- Divergence");
//       iter++;
//       break;
//     }

//     // Solve for the deltas - step 1 (ksi)
//     bk[0] = alpha + 2.*A[0][1];
//     ck[0] = -A[0][1];

//     fkx[0] = alpha*omega*L_phi_x[0][1];
//     fky[0] = alpha*omega*L_phi_y[0][1];
    
//     for(int i=2;i<TE-2;i++){
//       ak[i-2] = -A[0][i];
//       bk[i-1] = alpha + 2.*A[0][i];
//       ck[i-1] = -A[0][i];

//       fkx[i-1] = alpha*omega*L_phi_x[0][i];
//       fky[i-1] = alpha*omega*L_phi_y[0][i];
//     }

//     ak[TE-4] = -A[0][TE-2];
//     bk[TE-3] = alpha + 2.*A[0][TE-2];

//     fkx[TE-3] = alpha*omega*L_phi_x[0][TE-2];
//     fky[TE-3] = alpha*omega*L_phi_y[0][TE-2];
    
//     tridiagonal_matrix_solver(TE,ak,bk,ck,fkx,ukx);
//     tridiagonal_matrix_solver(TE,ak,bk,ck,fky,uky);

//     for(int i=1;i<TE-1;i++){
//       fx[0][i] = ukx[i-1];
//       fy[0][i] = uky[i-1];
//     }

//     for(int j=1;j<m-1;j++){
//       bk[0] = alpha + 2.*A[j][1];
//       ck[0] = -A[j][1];

//       fkx[0] = alpha*omega*L_phi_x[j][1];
//       fky[0] = alpha*omega*L_phi_y[j][1];
      
//       for(int i=2;i<n-2;i++){
//         ak[i-2] = -A[j][i];
//         bk[i-1] = alpha + 2.*A[j][i];
//         ck[i-1] = -A[j][i];

//         fkx[i-1] = alpha*omega*L_phi_x[j][i];
//         fky[i-1] = alpha*omega*L_phi_y[j][i];
//       }

//       ak[n-4] = -A[j][n-2];
//       bk[n-3] = alpha + 2.*A[j][n-2];

//       fkx[n-3] = alpha*omega*L_phi_x[j][n-2];
//       fky[n-3] = alpha*omega*L_phi_y[j][n-2];
      
//       tridiagonal_matrix_solver(n-2,ak,bk,ck,fkx,ukx);
//       tridiagonal_matrix_solver(n-2,ak,bk,ck,fky,uky);

//       for(int i=1;i<n-1;i++){
//         fx[j][i] = ukx[i-1];
//         fy[j][i] = uky[i-1];
//       }
//     }

//     print_2d_array_to_file(m,n,fx,"./mat_fx.dat",1);
//     print_2d_array_to_file(m,n,fy,"./mat_fy.dat",1);
//     exit(1);

//     // Solve for the deltas - step 2 (eta)
//     for(int i1=1,i2=n-2,k;i1<TE;i1++,i2--){
//       k = 0;
//       // printf("%d | %f %f\n",i,x[0][i],y[0][i]);
      
//       bn[k] = alpha + 2.*C[m-2][i2];
//       cn[k] = -C[m-2][i2];

//       fnx[k] = fx[m-2][i2];
//       fny[k] = fy[m-2][i2];
      
//       k++;

//       for(int j=m-3;j>=1;j--,k++){
//         an[k-1] = -C[j][i2];
//         bn[k] = alpha + 2.*C[j][i2];
//         cn[k] = -C[j][i2];

//         fnx[k] = fx[j][i2];
//         fny[k] = fy[j][i2];
//         // printf("%i\n",j);
//       }

//       for(int j=0;j<m-2;j++,k++){
//         an[k-1] = -C[j][i1];
//         bn[k] = alpha + 2.*C[j][i1];
//         cn[k] = -C[j][i1];

//         fnx[k] = fx[j][i1];
//         fny[k] = fy[j][i1];
//         // printf("%i\n",j);
//       }

//       an[m_wake-4] = -C[m-2][i1];
//       bn[m_wake-3] = alpha + 2.*C[m-2][i1];

//       fnx[m_wake-3] = fx[m-2][i1];
//       fny[m_wake-3] = fy[m-2][i1];

//       // print_1d_array(m-3,an);
//       // print_1d_array(m-2,bn);
//       // print_1d_array(m-3,cn);
//       // print_1d_array(m-2,fnx);
//       // print_1d_array(m-2,fny);

//       tridiagonal_matrix_solver(m_wake-2,an,bn,cn,fnx,unx);
//       tridiagonal_matrix_solver(m_wake-2,an,bn,cn,fny,uny);

//       for(int j=1,k1=(m_wake-1)/2,k2=k1-2;j<m-1;j++,k1++,k2--){
//         // printf("%02d %02d %02d\n",j,k1,k2);
//         Delta_x[j][i1] = unx[k1];
//         Delta_y[j][i1] = uny[k1];
//         Delta_x[j][i2] = unx[k2];
//         Delta_y[j][i2] = uny[k2];
//       }

//       // disp((m_wake-3)/2);
//       Delta_x[0][i1] = unx[(m_wake-1)/2];
//       Delta_y[0][i1] = uny[(m_wake-1)/2];
//       Delta_x[0][i1] = unx[(m_wake-1)/2];
//       Delta_y[0][i1] = uny[(m_wake-1)/2];
//     }

//     for(int i=TE;i<n-TE;i++){
//       bn[0] = alpha + 2.*C[1][i];
//       cn[0] = -C[1][i];

//       fnx[0] = fx[1][i];
//       fny[0] = fy[1][i];
      
//       for(int j=2;j<m-2;j++){
//         an[j-2] = -C[j][i];
//         bn[j-1] = alpha + 2.*C[j][i];
//         cn[j-1] = -C[j][i];

//         fnx[j-1] = fx[j][i];
//         fny[j-1] = fy[j][i];
//       }

//       an[m-4] = -C[m-2][i];
//       bn[m-3] = alpha + 2.*C[m-2][i];

//       fnx[m-3] = fx[m-2][i];
//       fny[m-3] = fy[m-2][i];

//       tridiagonal_matrix_solver(m-2,an,bn,cn,fnx,unx);
//       tridiagonal_matrix_solver(m-2,an,bn,cn,fny,uny);

//       for(int j=1;j<m-1;j++){
//         Delta_x[j][i] = unx[j-1];
//         Delta_y[j][i] = uny[j-1];
//       }
      
//     }

//      // Calculate new x and y
//      for(int j=1;j<m-1;j++){
//       for(int i=0;i<n-1;i++){
//         x[j][i] += Delta_x[j][i];
//         y[j][i] += Delta_y[j][i];
//       }
//     }

//     // Reapply "periodicity" boundary condition
//     for(int i1=1,i2=n-2;i1<TE;i1++,i2--){
//       x[0][i2] = x[0][i1];
//       y[0][i2] = y[0][i1];
//     }

//     if(!config->save_last_only){
//       save_results_qtimes(m,n,x,&iter,config,buffer,filename_save_x,
//                           &str_end_idx);
//       save_results_qtimes(m,n,y,&iter,config,buffer,filename_save_y,
//                           &str_end_idx);
//     }
//   }

//   // Save last iteration if it wasn't saved
//   if(iter%config->qtimes != 0 || config->save_last_only){
//     iter--; // To get the correct iteration number
//     sprintf(buffer,"L");
//     save_results_qtimes(m,n,x,&iter,config,buffer,filename_save_x,&str_end_idx);
//     sprintf(buffer,"L");
//     save_results_qtimes(m,n,y,&iter,config,buffer,filename_save_y,&str_end_idx);
//   }

//   free(L_phi_x);
//   free(L_phi_y);
//   free(Delta_x);
//   free(Delta_y);
//   free(A);
//   free(B);
//   free(C);
//   // free(D);
//   free(fx);
//   free(fy);
//   free(ak);
//   free(bk);
//   free(ck);
//   free(ukx);
//   free(uky);
//   free(fkx);
//   free(fky);
//   free(an);
//   free(bn);
//   free(cn);
//   free(unx);
//   free(uny);
//   free(fnx);
//   free(fny);
//   free(filename_save_x);
//   free(filename_save_y);
//   free(buffer);
//   free(filename_log_x);
//   free(filename_log_y);
//   fclose(file_log_x);
//   fclose(file_log_y);
// }

/*
kind of unnecessary? but i wanted to keep stuff standardized
- gigiaero, 20/05/2026, 1449 hours

okay, turns out i was incredibly wrong with the indexing here
- gigiaero, 21/05/2026, 1420 hours
*/
double uniform_scheme_der1_o2_central_ecm(int m,int n,double phi[m][n],int i){
  return (phi[1][i] - phi[1][n-i-1])*.5;
}

/*
- gigiaero, 20/05/2026, 1449 hours

indexing corrected, as noted above
- gigiaero, 21/05/2026, 1420 hours
*/
double uniform_scheme_der2_o2_central_ecm(int m,int n,double phi[m][n],int i,
                                          int axis){
  switch(axis){
    case 2: // Vertical
      return phi[1][i] - 2.*phi[0][i] + phi[1][n-i-1];
      break;

    case 3: // Mixed
      return (phi[1][i+1] - phi[1][n-i-1] - phi[1][i-1] + phi[1][n-i-1])*.25;
      break;

    default:
      puts("uniform_scheme_der2_o2_central: invalid axis");
      exit(15);
  }
}