#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/2d_arrays.h"
#include"../include/bi_air_lib.h"
#include"../include/eom_lib.h"
#include"../include/pom_lib.h"
#include"../include/strings.h"

#define div_ref 1e100
#define pi 3.14159265358979323846261003

/*
- gigiaero, 19/05/2026, 0951 hours
*/
void alpha_sequence(double *alpha,int *k,int iter,sim_prmtrs *config){
  if(config->alpha_seq){
    *alpha = config->alpha_H*pow(config->alpha/config->alpha_H,
             ((double)((*k) - 1))/((double)(config->M - 1)));

    (*k)++;

    if(*k > config->M)
      *k = 1;
  }
  else{
    if(iter == 1)
      *alpha = config->alpha;
  }
}

/*
- gigiaero, 20/05/2026, 0853 hours
*/
void alpha_sequence_aH(int m,int n,double x[m][n],double y[m][n],
                       sim_prmtrs *config,int k){
  if(k != config->M && config->set_alpha_H != 0 && config->alpha_seq){
    switch(config->set_alpha_H){
      case 1: // ADI
        config->alpha_H = min_physical_spacing(m,n,x,y);
        config->alpha_H = 4./config->alpha_H/config->alpha_H;
        break;
      
      case 2: // AF2
        config->alpha_H = 1./min_physical_spacing(m,n,x,y);
        break;

      default:
        puts("alpha_sequence_aH: invalid set_alpha_H");
        exit(11);
    }
  }
}

/*
- gigiaero, 27/04/2026, 1310 hours
*/
void calc_A(int m,int n,double A[m][n],double x[m][n],double y[m][n]){
  for(int j=1;j<m-1;j++){
    for(int i=0;i<n-1;i++)
      A[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,2),2.) + 
                pow(uniform_scheme_der1_o2_central(m,n,y,i,j,2),2.);
  }
}

/*
- gigiaero, 27/04/2026, 1310 hours

added correction for the periodicity condition
- gigiaero, 28/04/2026, 2105 hours
*/
void calc_B(int m,int n,double B[m][n],double x[m][n],double y[m][n]){
  for(int j=1;j<m-1;j++){
    B[j][0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j)*
              uniform_scheme_der1_o2_central(m,n,x,0,j,2) + 
              uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j)*
              uniform_scheme_der1_o2_central(m,n,y,0,j,2);
    // B[j][0] = (x[j][1] - x[j][n-2])*.5*
    //           uniform_scheme_der1_o2_central(m,n,x,0,j,2) + 
    //           (y[j][1] - y[j][n-2])*.5*
    //           uniform_scheme_der1_o2_central(m,n,y,0,j,2);

    for(int i=1;i<n-1;i++)
      B[j][i] = uniform_scheme_der1_o2_central(m,n,x,i,j,1)*
                uniform_scheme_der1_o2_central(m,n,x,i,j,2) + 
                uniform_scheme_der1_o2_central(m,n,y,i,j,1)*
                uniform_scheme_der1_o2_central(m,n,y,i,j,2);
  }
}

/*
- gigiaero, 27/04/2026, 1310 hours

added correction for the periodicity condition
- gigiaero, 28/04/2026, 2109 hours
*/
void calc_C(int m,int n,double C[m][n],double x[m][n],double y[m][n]){
  for(int j=1;j<m-1;j++){
    C[j][0] = pow(uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j),2.) + 
              pow(uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j),2.);
    // C[j][0] = pow((x[j][1] - x[j][n-2])*.5,2) + 
    //           pow((y[j][1] - y[j][n-2])*.5,2);

    for(int i=1;i<n-1;i++)
      C[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,1),2.) + 
                pow(uniform_scheme_der1_o2_central(m,n,y,i,j,1),2.);
  }
}

/*
- gigiaero, 27/04/2026, 1310 hours

added correction for the periodicity condition
- gigiaero, 28/04/2026, 2111 hours
*/
void calc_D(int m,int n,double D[m][n],double x[m][n],double y[m][n]){
  for(int j=1;j<m-1;j++){
    D[j][0] = pow(uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j)*
                  uniform_scheme_der1_o2_central(m,n,y,0,j,2) - 
                  uniform_scheme_der1_o2_central(m,n,x,0,j,2)*
                  uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j),2.);

    for(int i=1;i<n-1;i++)
      D[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,1)*
                    uniform_scheme_der1_o2_central(m,n,y,i,j,2) - 
                    uniform_scheme_der1_o2_central(m,n,x,i,j,2)*
                    uniform_scheme_der1_o2_central(m,n,y,i,j,1),2.);
  }
}

/*
- gigiaero, 07/05/2026, 2150 hours
*/
double control_P(control_prmtrs *ctrl,int ksi,int eta){
  double sum_l = 0.,sum_m = 0.;

  for(int l=0;l<ctrl->L;l++)
    sum_l += ctrl->al[l]*sgn(ksi - ctrl->ksi_l[l])*
                         exp(-ctrl->cl[l]*((double) abs(ksi - ctrl->ksi_l[l])));

  for(int m=0;m<ctrl->M;m++)
    sum_m += ctrl->bm[m]*sgn(ksi - ctrl->ksi_m[m])*
                         exp(-ctrl->dm[m]*sqrt(
                         pow((double) (ksi - ctrl->ksi_m[m]),2.) + 
                         pow((double) (eta - ctrl->eta_m[m]),2.)));

  return -sum_l - sum_m;
}

/*
- gigiaero, 07/05/2026, 2150 hours
*/
double control_Q(control_prmtrs *ctrl,int ksi,int eta){
  double sum_l = 0.,sum_m = 0.;

  for(int l=0;l<ctrl->L;l++)
    sum_l += ctrl->al[l]*sgn(eta - ctrl->eta_l[l])*
                         exp(-ctrl->cl[l]*((double) abs(eta - ctrl->eta_l[l])));

  for(int m=0;m<ctrl->M;m++)
    sum_m += ctrl->bm[m]*sgn(eta - ctrl->eta_m[m])*
                         exp(-ctrl->dm[m]*sqrt(
                         pow((double) (ksi - ctrl->ksi_m[m]),2.) + 
                         pow((double) (eta - ctrl->eta_m[m]),2.)));

  return -sum_l - sum_m;
}

/*
- gigiaero, 24/04/2026, 1352 hours

corrected a little mistake regarding th

- gigiaero, 29/04/2026, 2136 hours
*/
void cosspace(double *x,double xi,double xf,int n,int half){
  double midp;
  double th_i = pi/(((double) n) - 1.);
  double th;

  if(half){
    midp = xf - xi;
    th_i /= 2.;
  }
  else 
    midp = (xf - xi)/2.;
  
  th = th_i;

  x[0] = xi;
  for(int i=1;i<n;i++){
      x[i] = xi + midp*(1. - cos(th));
      th += th_i;
  }
}

void cst_airfoil(int n_pts,double *x,double *yu,double *yl,double *prmtrs,
                 double c){
  int n = (int) prmtrs[0];
  double *v_ex = malloc(sizeof(double)*(n+2));
  double *v_in = malloc(sizeof(double)*(n+2));
  double *A1 = malloc(sizeof(double)*(n+1));
  double *A2 = malloc(sizeof(double)*(n+1));
  double N1 = .5,N2 = 1.;
  double sum1,sum2,K;

  for(int i=0;i<n+2;i++){
    v_ex[i] = prmtrs[i+1];
    v_in[i] = prmtrs[i+1+n+2];
  }

  A1[0] = pow(2.*v_ex[0],.5);
  A2[0] = pow(2.*v_in[0],.5);
  for(int k=1;k<n;k++){
    A1[k] = v_ex[k];
    A2[k] = v_in[k];
  }
  A1[n] = tan(v_ex[n]) + v_ex[n+1]/c;
  A2[n] = tan(v_in[n]) + v_in[n+1]/c;

  for(int i=0;i<n_pts;i++){
    sum1 = 0.;
    sum2 = 0.;
    for(int r=0;r<n+1;r++){
      // K is just the binomial coefficient, but bear with me here
      K = factorial(n)/(factorial(r)*factorial(n-r))*
          pow(x[i],r)*pow(1. - x[i],n - r);
      sum1 += A1[r]*K;
      sum2 += A2[r]*K;
    }

    yu[i] = pow(x[i],N1)*pow(1.-x[i],N2)*sum1;
    yl[i] = pow(x[i],N1)*pow(1.-x[i],N2)*sum2;
  }

  free(v_ex);
  free(v_in);
  free(A1);
  free(A2);
}

/*
- gigiaero, 25/04/2026, 1134 hours

function inputs updated
- gigiaero, 29/05/2026, 1336 hours
*/
void cst_prmtrs(msh_prmtrs *msh){
  int n = (int) msh->af_prmtrs[0];
  double *v_ex = malloc(sizeof(double)*(n+2));
  double *v_in = malloc(sizeof(double)*(n+2));
  char *filename = malloc(sizeof(char)*100);
  sprintf(filename,"../P02/reverse_cst/");
  FILE *input;
  
  switch(msh->cst_foil){
    case 1: // n = 10
      strcat(filename,"whitcomb_cst_prmtrs.dat");
      break;
    
    case 2: // n = 10
      strcat(filename,"fx61-163_cst_prmtrs.dat");
      break;
    
    case 3: // n = 10
      strcat(filename,"naca64A005.92_cst_prmtrs.dat");
      break;
    
    case 4: // n = 10
      strcat(filename,"s1223_cst_prmtrs.dat");
      break;

    default:
      puts("cst_prmtrs: invalid airfoil")  ;
      exit(108);
  }

  input = fopen(filename,"r");

  for(int i=0;i<n+2;i++){
    if(fscanf(input,"%lf",&v_ex[i]) == EOF)
      break;
  }

  for(int i=0;i<n+2;i++){
    if(fscanf(input,"%lf",&v_in[i]) == EOF)
      break;
  }

  for(int i=0;i<n+2;i++){
    msh->af_prmtrs[i+1] = v_ex[i];
    msh->af_prmtrs[i+1+n+2] = v_in[i];
  }

  free(v_ex);
  free(v_in);
  free(filename);
  fclose(input);
}

/*
prmtrs = {rx,ry,dx,dy}
also makes circles when rx = ry
- gigiaero, 24/04/2026,1635 hours
*/
void ellipse(double *x,double *y,double *prmtrs,int n,int invert_th){
  double *th = malloc(sizeof(double)*n);

  if(invert_th)
    linspace(th,2.*pi,0.,n);
  else
    linspace(th,0.,2.*pi,n);

  for(int i=0;i<n;i++){
    x[i] = prmtrs[2] + prmtrs[0]*cos(th[i]);
    y[i] = prmtrs[3] + prmtrs[1]*sin(th[i]);
  }
    
  free(th);
}

/*
- gigiaero, 26/04/2026, 1003 hours
*/
void evaluate_delta_form_eom(sim_prmtrs *config,msh_prmtrs *msh,
                             control_prmtrs *c_prmtrs,int init_only){
  double (*x)[msh->IMAX] = calloc(msh->JMAX,sizeof *x);
  double (*y)[msh->IMAX] = calloc(msh->JMAX,sizeof *y);

  save_prmtrs_sim(config);
  save_prmtrs_msh(config->casename,msh);

  initialize_mesh(msh->JMAX,msh->IMAX,x,y,msh);

  if(config->save_i_c){
    char *filename = malloc(sizeof(char)*200);

    sprintf(filename,"%s%s",config->casename,"_x_initial.dat");
    print_2d_array_to_file(msh->JMAX,msh->IMAX,x,filename,0);
    sprintf(filename,"%s%s",config->casename,"_y_initial.dat");
    print_2d_array_to_file(msh->JMAX,msh->IMAX,y,filename,0);

    free(filename);
  }

  if(!init_only){
    switch(config->Ntype){
      case 1:
        puts("SLOR, elliptical O mesh");
        solve_slor_2d_rectangular_eom(msh->JMAX,msh->IMAX,x,y,config,c_prmtrs);
        break;
      
      case 2:
        puts("ADI, elliptical O mesh");
        solve_adi_2d_rectangular_eom(msh->JMAX,msh->IMAX,x,y,config,c_prmtrs);
        break;
      
      case 3:
        puts("AF2, elliptical O mesh");
        solve_af2_2d_rectangular_eom(msh->JMAX,msh->IMAX,x,y,config,c_prmtrs);
        break;

      case 4:
        puts("ADI, elliptical O mesh, nonperiodic");
        solve_adi_2d_rectangular_eom_np(msh->JMAX,msh->IMAX,x,y,config,
                                        c_prmtrs);
        break;
        
      default:
        puts("evaluate_delta_form_eom: Invalid Ntype");
        exit(32);
    }
  }

  free(x);
  free(y);
}

// /*
// for usage in the full potential code
// - gigiaero, 28/05/2026, 2057 hours
// */
// void evaluate_delta_form_eom_fp(int m,int n,double x[m][n],double y[m][n],
//                                 sim_prmtrs *config,msh_prmtrs *msh,
//                                 control_prmtrs *c_prmtrs,int init_only){
//   initialize_mesh(msh->JMAX,msh->IMAX,x,y,msh);

//   if(init_only){
//     char *filename = malloc(sizeof(char)*200);
//     sprintf(filename,"%s%s",config->casename,"_x_initial.dat");
//     print_2d_array_to_file(msh->JMAX,msh->IMAX,x,filename,0);
//     sprintf(filename,"%s%s",config->casename,"_y_initial.dat");
//     print_2d_array_to_file(msh->JMAX,msh->IMAX,y,filename,0);
//     free(filename);
//   }
//   else{
//     config->save_last_only = 1;
//     config->qtimes = 9999999999;

//     switch(config->Ntype){
//       case 1:
//         puts("SLOR, elliptical O mesh");
//         solve_slor_2d_rectangular_eom(msh->JMAX,msh->IMAX,x,y,config,c_prmtrs);
//         break;
      
//       case 2:
//         puts("ADI, elliptical O mesh");
//         solve_adi_2d_rectangular_eom(msh->JMAX,msh->IMAX,x,y,config,c_prmtrs);
//         break;
      
//       case 3:
//         puts("AF2, elliptical O mesh");
//         solve_af2_2d_rectangular_eom(msh->JMAX,msh->IMAX,x,y,config,c_prmtrs);
//         break;

//       case 4:
//         puts("ADI, elliptical O mesh, nonperiodic");
//         solve_adi_2d_rectangular_eom_np(msh->JMAX,msh->IMAX,x,y,config,
//                                         c_prmtrs);
//         break;
        
//       default:
//         puts("evaluate_delta_form_eom: Invalid Ntype");
//         exit(32);
//     }
//   }
// }

/*
- gigiaero, 24/04/2026, 2231 hours
*/
void init_af_bi_air(double *x,double *y,double *x_axis,int chord_n,
                    msh_prmtrs *msh){
  x[chord_n] = x_axis[0];
  y[chord_n] = x_axis[0];
  for(int i1=chord_n-2,i2=chord_n,k=1;i1>=0;i1--,i2++,k++){
    x[i1] = x_axis[k];
    x[i2] = x_axis[k];
    y[i1] = -bi_air_shape(msh->af_prmtrs[0],x_axis[k]);
    y[i2] = bi_air_shape(msh->af_prmtrs[0],x_axis[k]);
  }
}

void init_af_cst(double *x,double *y,double *x_axis,int chord_n,
                 msh_prmtrs *msh){
  double *yu = malloc(sizeof(double)*chord_n);
  double *yl = malloc(sizeof(double)*chord_n);
  
  cst_airfoil(chord_n,x_axis,yu,yl,msh->af_prmtrs,msh->c);
  
  y[chord_n] = yu[0];
  for(int i1=chord_n-2,i2=chord_n,k=1;i1>=0;i1--,i2++,k++){
    x[i1] = x_axis[k];
    x[i2] = x_axis[k];
    y[i1] = -yl[k];
    y[i2] = yu[k];
  }

  free(yu);
  free(yl);
}

void init_af_naca4(double *x,double *y,double *x_axis,int chord_n,
                   msh_prmtrs *msh){
  double *xu = malloc(sizeof(double)*chord_n);
  double *xl = malloc(sizeof(double)*chord_n);
  double *yu = malloc(sizeof(double)*chord_n);
  double *yl = malloc(sizeof(double)*chord_n);
  
  naca4(chord_n,x_axis,xu,xl,yu,yl,msh->af_prmtrs);
  
  x[chord_n] = xu[0];
  y[chord_n] = yu[0];
  for(int i1=chord_n-2,i2=chord_n,k=1;i1>=0;i1--,i2++,k++){
    x[i1] = xl[k];
    x[i2] = xu[k];
    y[i1] = yl[k];
    y[i2] = yu[k];
  }

  free(xu);
  free(xl);
  free(yu);
  free(yl);
}

/*
- gigiaero, 24/04/2026, 2241 hours
*/
void init_type1(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  double *vrx = malloc(sizeof(double)*msh->JMAX);
  double *vry = malloc(sizeof(double)*msh->JMAX);

  linspace(vrx,msh->c*.5,msh->end_prmtrs[0],msh->JMAX);
  linspace(vry,max_thickness(m,n,y)*.5,msh->end_prmtrs[1],msh->JMAX);

  for(int j=1;j<msh->JMAX-1;j++){
    msh->end_prmtrs[0] = vrx[j];
    msh->end_prmtrs[1] = vry[j];
    ellipse(x[j],y[j],msh->end_prmtrs,msh->IMAX,1);
  }

  // restore original values
  msh->end_prmtrs[0] = vrx[msh->JMAX-1];
  msh->end_prmtrs[1] = vry[msh->JMAX-1];

  free(vrx);
  free(vry);
}

/*
- gigiaero, 24/04/2026, 2204 hours

changed i loop as done for init_type3
- gigiaero, 21/05/2026, 1458 hours
*/
void init_type2(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  double *vx = malloc(sizeof(double)*msh->JMAX);
  double *vy = malloc(sizeof(double)*msh->JMAX);

  for(int i=0;i<n;i++){
    linspace(vx,x[0][i],x[msh->JMAX-1][i],msh->JMAX);
    linspace(vy,y[0][i],y[msh->JMAX-1][i],msh->JMAX);
    
    for(int j=1;j<msh->JMAX-1;j++){
      x[j][i] = vx[j];
      y[j][i] = vy[j];
    }
  }

  free(vx);
  free(vy);
}

/*
- gigiaero, 24/04/2026, 2204 hours

in the i loop, the bound was changed from msh->IMAX to n to accomodate the ecm
solution
- gigiaero, 20/05/2026, 1412 hours
*/
void init_type3(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  double *vx = malloc(sizeof(double)*msh->JMAX);
  double *vy = malloc(sizeof(double)*msh->JMAX);

  for(int i=0;i<n;i++){
    cosspace(vx,x[0][i],x[msh->JMAX-1][i],msh->JMAX,1);
    cosspace(vy,y[0][i],y[msh->JMAX-1][i],msh->JMAX,1);
    
    for(int j=1;j<msh->JMAX-1;j++){
      x[j][i] = vx[j];
      y[j][i] = vy[j];
    }
  }

  free(vx);
  free(vy);
}

// // parabolic mesh
// void init_type4(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  
// }

/*
- gigiaero, 24/04/2026, 2204 hours
*/
void initialize_mesh(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  int chord_n = (n+1)/2;
  double *x_axis = malloc(sizeof(double)*chord_n);
  cosspace(x_axis,0.,msh->c,chord_n,0);
  
  switch(msh->af_type){
    case 1:
      init_af_bi_air(x[0],y[0],x_axis,chord_n,msh);
      break;

    case 2:
      init_af_naca4(x[0],y[0],x_axis,chord_n,msh);
      break;

    case 3:
      init_af_cst(x[0],y[0],x_axis,chord_n,msh);
      break;

    default:
      puts("initialize_mesh: invalid af_type");
      exit(13);
  }

  ellipse(x[m-1],y[m-1],msh->end_prmtrs,n,1);

  switch(msh->init_type){
    case 1:
      init_type1(m,n,x,y,msh);
      break;

    case 2:
      init_type2(m,n,x,y,msh);
      break;

    case 3:
      init_type3(m,n,x,y,msh);
      break;

    case 4:
      parabolic_mesh(m,n,x,y,msh);
      break;

    default:
      puts("initialize_mesh: invalid init_type");
      exit(4);
  }

  free(x_axis);
}

/*
- gigiaero, 27/04/2026, 1317 hours

added correction for the periodicity condition
- gigiaero, 28/04/2026, 2115 hours

added control functions P and Q
- gigiaero, 07/05/2026, 2150 hours
*/
void L_phi_eom(int m,int n,double L_phi_x[m][n],double L_phi_y[m][n],
               double x[m][n],double y[m][n],double A[m][n],
               double B[m][n],double C[m][n],double D[m][n],
               control_prmtrs *c_prmtrs){
  double P,Q;

  for(int j=1;j<m-1;j++){
    P = control_P(c_prmtrs,0,j);
    Q = control_Q(c_prmtrs,0,j);

    L_phi_x[j][0] =A[j][0]*uniform_scheme_der2_o2_central_prdc_ksi(m,n,x,j,1) - 
                   2.*B[j][0]*uniform_scheme_der2_o2_central_prdc_ksi(m,n,x,j,3)
                   +C[j][0]*uniform_scheme_der2_o2_central(m,n,x,0,j,2) + 
                    D[j][0]*(P*uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j)
                           + Q*uniform_scheme_der1_o2_central(m,n,x,0,j,2));
                    
    L_phi_y[j][0] =A[j][0]*uniform_scheme_der2_o2_central_prdc_ksi(m,n,y,j,1) - 
                   2.*B[j][0]*uniform_scheme_der2_o2_central_prdc_ksi(m,n,y,j,3)
                   +C[j][0]*uniform_scheme_der2_o2_central(m,n,y,0,j,2) +
                    D[j][0]*(P*uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j)
                           + Q*uniform_scheme_der1_o2_central(m,n,y,0,j,2));

    for(int i=1;i<n-1;i++){
      P = control_P(c_prmtrs,i,j);
      Q = control_Q(c_prmtrs,i,j);
      L_phi_x[j][i] = A[j][i]*uniform_scheme_der2_o2_central(m,n,x,i,j,1) - 
                       2.*B[j][i]*uniform_scheme_der2_o2_central(m,n,x,i,j,3) + 
                       C[j][i]*uniform_scheme_der2_o2_central(m,n,x,i,j,2) +
                       D[j][i]*(P*uniform_scheme_der1_o2_central(m,n,x,i,j,1) +
                                Q*uniform_scheme_der1_o2_central(m,n,x,i,j,2));
                       
      L_phi_y[j][i] = A[j][i]*uniform_scheme_der2_o2_central(m,n,y,i,j,1) - 
                       2.*B[j][i]*uniform_scheme_der2_o2_central(m,n,y,i,j,3) + 
                       C[j][i]*uniform_scheme_der2_o2_central(m,n,y,i,j,2) + 
                       D[j][i]*(P*uniform_scheme_der1_o2_central(m,n,y,i,j,1) +
                                Q*uniform_scheme_der1_o2_central(m,n,y,i,j,2));
    }

    L_phi_x[j][n-1] = L_phi_x[j][0];
    L_phi_y[j][n-1] = L_phi_y[j][0];
  }
}

/*
- gigiaero, 24/04/2026, 1313 hours
*/
void linspace(double *x,double xi,double xf,int n){

  x[0] = xi;
  double step = (xf - xi)/((double) n - 1);

  for(int i=1;i<n;i++)
    x[i] = x[i-1] + step;
}

/*
- gigiaero, 07/05/2026, 2112 hours
*/
void malloc_c_prmtrs(control_prmtrs *c_prmtrs){
  c_prmtrs->al = malloc(sizeof(double)*c_prmtrs->L);
  c_prmtrs->bm = malloc(sizeof(double)*c_prmtrs->M);
  c_prmtrs->cl = malloc(sizeof(double)*c_prmtrs->L);
  c_prmtrs->dm = malloc(sizeof(double)*c_prmtrs->M);
  c_prmtrs->ksi_l = malloc(sizeof(int)*c_prmtrs->L);
  c_prmtrs->ksi_m = malloc(sizeof(int)*c_prmtrs->M);
  c_prmtrs->eta_l = malloc(sizeof(int)*c_prmtrs->L);
  c_prmtrs->eta_m = malloc(sizeof(int)*c_prmtrs->M);
}

/*
- gigiaero, 20/05/2026, 0930 hours
*/
double min_physical_spacing(int m,int n,double x[m][n],double y[m][n]){
  // double kxx=0.,kxy=0.,kyx=0.,kyy=0.;
  double min=1e10;
  double deltaxy,deltaxx,deltayx,deltayy;

  for(int j=0;j<m-1;j++){
    for(int i=0;i<n-1;i++){
      deltaxx = fabs(x[j][i+1] - x[j][i]);
      deltaxy = fabs(x[j+1][i] - x[j][i]);
      deltayx = fabs(y[j][i+1] - y[j][i]);
      deltayy = fabs(y[j+1][i] - y[j][i]);

      /* limpar isto quando eu estiver certa que está tudo em ordem */
      // disp(deltaxx);
      // disp(deltaxy);
      // disp(deltayx);
      // disp(deltayy);

      if(deltaxx > deltaxy && deltaxy != 0)
        deltaxx = deltaxy;

      if(deltayy > deltayx && deltayx != 0)
        deltayy = deltayx;

      if(min > deltaxx && deltaxx != 0)
        min = deltaxx;
      
      if(min > deltayy && deltayy != 0)
        min = deltayy;
      
      // disp(min);
      // putchar('\n');
      // min = (deltaxx > deltayy) ? deltaxx : deltayy;
    }
  }

  // puts("-------------------------------"); putchar('\n');

  for(int j=0;j<m;j++){
    deltaxx = fabs(x[j][n-1] - x[j][n-2]);
    // deltaxy = x[j+1][i] - x[j][i];
    deltayx = fabs(y[j][n-1] - y[j][n-2]);
    // deltayy = y[j+1][i] - y[j][i];

    // disp(deltaxx);
    // disp(deltayx);

    if(min > deltaxy && deltaxy != 0)
      min = deltaxy;

    if(min > deltayx && deltayx != 0)
      min = deltayx;

    // min = (deltaxx > deltayx) ? deltaxx : deltayx; 
    // disp(min);
    // putchar('\n');
  }

  // puts("-------------------------------"); putchar('\n');

  for(int i=0;i<n;i++){
    // deltaxx = x[j][i+1] - x[j][i];
    deltaxy = fabs(x[m-1][i] - x[m-2][i]);
    // deltayx = y[j][i+1] - y[j][i];
    deltayy = fabs(y[m-1][i] - y[m-2][i]);

    // disp(deltaxy);
    // disp(deltayy);

    if(min > deltaxy)
      min = deltaxy;

    if(min > deltayx)
      min = deltayx;

    // min = (deltaxy > deltayy) ? deltaxy : deltayy;
    // disp(min); putchar('\n');
  }

  return min;
}

/*
- gigiaero, 24/04/2026, 2241 hours
*/
double max_thickness(int m,int n,double y[m][n]){
  int chord_n = (n+1)/2;
  double t = 0.,s;

  for(int i1=chord_n-2,i2=chord_n;i1>=0;i1--,i2++){
    s = fabs(y[0][i1] - y[0][i2]); 
    if(s > t)
      t = s;
  }

  return t;
}

/*
- gigiaero, 10/11/2023, 1526 hours
*/
void naca4(int n,double *x,double *xu,double *xl,double *yu,double *yl,
           double *prmtrs){
  double *thcknss = malloc(sizeof(double)*n);

  double m = prmtrs[0]/100.;
  double p = prmtrs[1]/10.;
  double t = prmtrs[2]/100.;

  double a0 = 0.2969;
  double a1 = -0.1260;
  double a2 = -0.3516;
  double a3 = 0.2843;
  // double a4 = -0.1015; // open trailing edge
  double a4 = -0.1036;

  for(int i=0;i<n;i++)
    thcknss[i] = 5.*t*(a0*pow(x[i],0.5) + a1*x[i] + a2*pow(x[i],2.) +
                       a3*pow(x[i],3.) + a4*pow(x[i],4.));

  // Symmetric 
  if(m==0. && p==0.){
    for(int i=0;i<n;i++){
      xu[i] = x[i];
      yu[i] = thcknss[i];
      xl[i] = x[i];
      yl[i] = -thcknss[i];
    }
  }
  // Asymmetric
  else{
    double *crvtr = malloc(sizeof(double)*n);
    double *slope = malloc(sizeof(double)*n);
    double *theta = malloc(sizeof(double)*n);

    for(int i=0;i<n;i++){
      if(x[i] >= 0. && x[i] <= p)
        crvtr[i] = m/pow(p,2.)*(2.*p*x[i] - pow(x[i],2.));

      else if(x[i] >= p && x[i] <= 1)
        crvtr[i] = m/pow(1.-p,2.)*((1. - 2.*p) + 2.*p*x[i] - pow(x[i],2.));
    }

    for(int i=0;i<n;i++){
      if(x[i] >= 0. && x[i] <= p)
        slope[i] = 2.*m/pow(p,2.)*(p - x[i]);

      else if(x[i] >= p && x[i] <= 1.)
        slope[i] = 2.*m/pow(1.-p,2.)*(p - x[i]);

      theta[i] = atan(slope[i]);
    }

    for(int i=0;i<n;i++){
      xu[i] = x[i] - thcknss[i]*sin(theta[i]);
      yu[i] = crvtr[i] + thcknss[i]*cos(theta[i]);
      xl[i] = x[i] + thcknss[i]*sin(theta[i]);
      yl[i] = crvtr[i] - thcknss[i]*cos(theta[i]);
    }

    free(crvtr);
    free(slope);
    free(theta);
  }
	
  free(thcknss);
}

/*
- gigiaero, 29/05/2026, 1342 hours
*/
void save_prmtrs_msh(char* casename,msh_prmtrs *msh){
  char *filename = malloc(sizeof(char)*200);
  FILE *output;

  sprintf(filename,"%s_prmtrs_msh.dat",casename);

  output = fopen(filename,"w");

  fprintf(output,"/* IMAX */\n");
  fprintf(output,"msh.IMAX = %d\n",msh->IMAX);
  fprintf(output,"/* JMAX */\n");
  fprintf(output,"msh.JMAX = %d\n",msh->JMAX);
  fprintf(output,"/* c */\n");
  fprintf(output,"msh.c = %f\n",msh->c);
  fprintf(output,"/* end_prmtrs */\n");
  fprintf(output,"msh.end_prmtrs[0] = %f\n",msh->end_prmtrs[0]);
  fprintf(output,"msh.end_prmtrs[1] = %f\n",msh->end_prmtrs[1]);
  fprintf(output,"msh.end_prmtrs[2] = %f\n",msh->end_prmtrs[2]);
  fprintf(output,"msh.end_prmtrs[3] = %f\n",msh->end_prmtrs[3]);
  fprintf(output,"/* init_type */\n");
  fprintf(output,"msh.init_type = %d\n",msh->init_type);
  fprintf(output,"/* af_type */\n");
  fprintf(output,"msh.af_type = %d\n",msh->af_type);
  if(msh->af_type == 1){
    fprintf(output,"/* af_prmtrs (bi_air) */\n");
    fprintf(output,"msh.af_prmtrs[0] = %f\n",msh->af_prmtrs[0]);
  }
  else if(msh->af_type == 2){
    fprintf(output,"/* af_prmtrs (naca4) */\n");
    fprintf(output,"msh.af_prmtrs[0] = %f\n",msh->af_prmtrs[0]);
    fprintf(output,"msh.af_prmtrs[1] = %f\n",msh->af_prmtrs[1]);
    fprintf(output,"msh.af_prmtrs[2] = %f\n",msh->af_prmtrs[2]);
  }
  else if(msh->af_type == 3){
    fprintf(output,"/* af_prmtrs (cst) */\n");
    fprintf(output,"int n = %d\n",(int) msh->af_prmtrs[0]);
    fprintf(output,"int cst_foil = %d\n",msh->cst_foil);
  }
  else{
    puts("save_prmtrs_msh: invalid af_type");
    fclose(output);
    free(filename);
    exit(108);
  }
  
  fclose(output);
  free(filename);
}

/*
defined in the context of fullp
- gigiaero, 02/06/2026, 0920 hours
*/
double scheme_der1_o2_central_prdc_ksi(int m,int n,double phi[m][n],
                                       double x[m][n],int j){
  return (phi[j][1] - phi[j][n-2])/(x[j][1] - x[j][n-2]);
}

/*
- gigiaero, 22/05/2026, 1101
*/
void set_control_prmtrs(int c_type,control_prmtrs *c_prmtrs,msh_prmtrs *msh){
  switch(c_type){
    case 0: // deactivated
      c_prmtrs->L = 1;
      c_prmtrs->M = 1;
      malloc_c_prmtrs(c_prmtrs);
      c_prmtrs->al[0] = 0;
      c_prmtrs->bm[0] = 0;
      c_prmtrs->cl[0] = 1;
      c_prmtrs->dm[0] = 1;
      c_prmtrs->ksi_l[0] = 0;
      c_prmtrs->eta_l[0] = 0;
      c_prmtrs->ksi_m[0] = 0;
      c_prmtrs->eta_m[0] = 0;
      break;

    /*
    group near the airfoil and the wake
    whitcomb airfoil
    config.alpha_seq = 1;
    config.alpha = 1e-2;
    config.alpha_H = 1e3;
    config.set_alpha_H = 0;
    config.M = 10;
    */
    case 1: 
      c_prmtrs->L = 2;
      c_prmtrs->M = 1;
      malloc_c_prmtrs(c_prmtrs);
      c_prmtrs->al[0] = 20; c_prmtrs->al[1] = c_prmtrs->al[0];
      c_prmtrs->bm[0] = 0;
      c_prmtrs->cl[0] = .5; c_prmtrs->cl[1] = c_prmtrs->cl[0];
      c_prmtrs->dm[0] = 1;
      c_prmtrs->ksi_l[0] = 0; c_prmtrs->ksi_l[1] = msh->IMAX-1;
      c_prmtrs->eta_l[0] = 0; c_prmtrs->eta_l[1] = 0;
      c_prmtrs->ksi_m[0] = 0;
      c_prmtrs->eta_m[0] = 0;
      break;
    
    /*
    refine leading edge
    */
    case 2:
      c_prmtrs->L = 1;
      c_prmtrs->M = 1;
      malloc_c_prmtrs(c_prmtrs);
      c_prmtrs->al[0] = 0;
      c_prmtrs->bm[0] = 1000;
      c_prmtrs->cl[0] = 1;
      c_prmtrs->dm[0] = .7;
      c_prmtrs->ksi_l[0] = 0;
      c_prmtrs->eta_l[0] = 0;
      c_prmtrs->ksi_m[0] = (msh->IMAX-1)/2;
      c_prmtrs->eta_m[0] = 0;
      break;

    /*
    refine trailing edge
    */
    case 3:
      c_prmtrs->L = 1;
      c_prmtrs->M = 2;
      malloc_c_prmtrs(c_prmtrs);
      c_prmtrs->al[0] = 0;
      c_prmtrs->bm[0] = 1000; c_prmtrs->bm[1] = c_prmtrs->bm[0];
      c_prmtrs->cl[0] = 1;
      c_prmtrs->dm[0] = .7; c_prmtrs->dm[1] = c_prmtrs->dm[0];
      c_prmtrs->ksi_l[0] = 0;
      c_prmtrs->eta_l[0] = 0;
      c_prmtrs->ksi_m[0] = 0; c_prmtrs->ksi_m[1] = msh->IMAX-1;
      c_prmtrs->eta_m[0] = 0; c_prmtrs->eta_m[1] = c_prmtrs->eta_m[0]; 
      break;

    /*
    refine shockwave region and close to the airfoil
    */
    case 4:
      c_prmtrs->L = 2;
      c_prmtrs->M = 1;
      malloc_c_prmtrs(c_prmtrs);
      c_prmtrs->al[0] = 15; c_prmtrs->al[1] = c_prmtrs->al[0];
      c_prmtrs->bm[0] = 0;
      c_prmtrs->cl[0] = 1; c_prmtrs->cl[1] = c_prmtrs->cl[0]; 
      c_prmtrs->dm[0] = 1;
      c_prmtrs->ksi_l[0] = 20; c_prmtrs->ksi_l[1] = 93 - c_prmtrs->ksi_l[0] - 1;
      c_prmtrs->eta_l[0] = 0; c_prmtrs->eta_l[1] = 0;
      c_prmtrs->ksi_m[0] = 0;
      c_prmtrs->eta_m[0] = 0;
      break;

    /*
    refine underside near trailing edge (imax = 93)
    */
    case 5:
      c_prmtrs->L = 1;
      c_prmtrs->M = 1;
      malloc_c_prmtrs(c_prmtrs);
      c_prmtrs->al[0] = 0; 
      c_prmtrs->bm[0] = 200;
      c_prmtrs->cl[0] = 1; 
      c_prmtrs->dm[0] = .8;
      c_prmtrs->ksi_l[0] = 0; 
      c_prmtrs->eta_l[0] = 0; 
      c_prmtrs->ksi_m[0] = 12;
      c_prmtrs->eta_m[0] = 0;
      break;

    /*
    do everything at the same time (imax=93)
    */
    case 6:
      c_prmtrs->L = 4;
      c_prmtrs->M = 4;
      malloc_c_prmtrs(c_prmtrs);
      // refine near airfoil surface and wake
      c_prmtrs->al[0] = 15; c_prmtrs->al[1] = c_prmtrs->al[0];
      c_prmtrs->cl[0] = .8; c_prmtrs->cl[1] = c_prmtrs->cl[0];
      c_prmtrs->ksi_l[0] = 0; c_prmtrs->ksi_l[1] = msh->IMAX-1;
      c_prmtrs->eta_l[0] = 0; c_prmtrs->eta_l[1] = 0;
      // refine leading edge
      c_prmtrs->bm[0] = 1000;
      c_prmtrs->dm[0] = .7;
      c_prmtrs->ksi_m[0] = (msh->IMAX-1)/2;
      c_prmtrs->eta_m[0] = 0;
      // refine trailing edge
      c_prmtrs->bm[1] = 1000; c_prmtrs->bm[2] = c_prmtrs->bm[1];
      c_prmtrs->dm[1] = .7; c_prmtrs->dm[2] = c_prmtrs->dm[1];
      c_prmtrs->ksi_m[1] = 0; c_prmtrs->ksi_m[2] = msh->IMAX-1;
      c_prmtrs->eta_m[1] = 0; c_prmtrs->eta_m[2] = c_prmtrs->eta_m[1]; 
      // refine shockwave regions and close to the airfoil
      c_prmtrs->al[2] = 15; c_prmtrs->al[3] = c_prmtrs->al[2];
      c_prmtrs->cl[2] = 1; c_prmtrs->cl[3] = c_prmtrs->cl[2]; 
      c_prmtrs->ksi_l[2] = 20; c_prmtrs->ksi_l[3] = 93 - c_prmtrs->ksi_l[2] - 1;
      c_prmtrs->eta_l[2] = 0; c_prmtrs->eta_l[3] = 0;
      // refine underside near trailing edge
      c_prmtrs->bm[0] = 200;
      c_prmtrs->dm[0] = .8;
      c_prmtrs->ksi_m[0] = 12;
      c_prmtrs->eta_m[0] = 0;
      break;
    
    case 50: // ecm: refine wake near trailing edge
      c_prmtrs->L = 1;
      c_prmtrs->M = 6;
      malloc_c_prmtrs(c_prmtrs);
      c_prmtrs->al[0] = 0;
      c_prmtrs->bm[0] = 160; c_prmtrs->bm[1] = 160; 
      c_prmtrs->bm[2] = 160; c_prmtrs->bm[3] = 160;
      c_prmtrs->bm[4] = 160; c_prmtrs->bm[5] = 160;
      c_prmtrs->cl[0] = .8;
      c_prmtrs->dm[0] = .5; c_prmtrs->dm[1] = .5;
      c_prmtrs->dm[2] = .5; c_prmtrs->dm[3] = .5;
      c_prmtrs->dm[4] = .5; c_prmtrs->dm[5] = .5;
      c_prmtrs->ksi_l[0] = 0;
      c_prmtrs->eta_l[0] = 0;
      c_prmtrs->ksi_m[0] = msh->JMAX; 
      c_prmtrs->ksi_m[1] = msh->IMAX - msh->JMAX;
      c_prmtrs->ksi_m[2] = msh->JMAX-1; 
      c_prmtrs->ksi_m[3] = msh->IMAX - msh->JMAX+1;
      c_prmtrs->ksi_m[4] = msh->JMAX-2; 
      c_prmtrs->ksi_m[5] = msh->IMAX - msh->JMAX+2;
      c_prmtrs->eta_m[0] = 0; 
      c_prmtrs->eta_m[1] = 0;
      c_prmtrs->eta_m[2] = 0;
      c_prmtrs->eta_m[3] = 0;
      c_prmtrs->eta_m[4] = 0;
      c_prmtrs->eta_m[5] = 0;
      break;
    

    default:
      puts("set_control_prmtrs: invalid c_type");
      exit(118);
  }
}

/*
- gigiaero, 07/05/2026, 1751 hours
*/
void solve_adi_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config,control_prmtrs *c_prmtrs){
  // Solver variables
  double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
  double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
  double (*Delta_x)[n] = calloc(m,sizeof *Delta_x);
  double (*Delta_y)[n] = calloc(m,sizeof *Delta_y);
  double (*A)[n] = calloc(m,sizeof *A);
  double (*B)[n] = calloc(m,sizeof *B);
  double (*C)[n] = calloc(m,sizeof *C);
  double (*D)[n] = calloc(m,sizeof *D);
  double (*fx)[n] = calloc(m,sizeof *fx);
  double (*fy)[n] = calloc(m,sizeof *fy);
  double *ak = malloc(sizeof(double)*(n-1));
  double *bk = malloc(sizeof(double)*(n-1));
  double *ck = malloc(sizeof(double)*(n-1));
  double *fkx = malloc(sizeof(double)*(n-1));
  double *fky = malloc(sizeof(double)*(n-1));
  double *ukx = malloc(sizeof(double)*(n-1));
  double *uky = malloc(sizeof(double)*(n-1));
  double *an = malloc(sizeof(double)*(m-3));
  double *bn = malloc(sizeof(double)*(m-2));
  double *cn = malloc(sizeof(double)*(m-3));
  double *fnx = malloc(sizeof(double)*(m-2));
  double *fny = malloc(sizeof(double)*(m-2));
  double *unx = malloc(sizeof(double)*(m-2));
  double *uny = malloc(sizeof(double)*(m-2));
  double omega = config->w;
  double alpha;
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
    alpha_sequence_aH(m,n,x,y,config,k);
    alpha_sequence(&alpha,&k,iter,config);

    calc_A(m,n,A,x,y);
    calc_B(m,n,B,x,y);
    calc_C(m,n,C,x,y);
    calc_D(m,n,D,x,y);

    // Calculate residual operators
    L_phi_eom(m,n,L_phi_x,L_phi_y,x,y,A,B,C,D,c_prmtrs);

    res_x = 0.;
    res_y = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
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

    // Solve for the deltas - step 1 (ksi)
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        ak[i] = -A[j][i];
        bk[i] = alpha + 2.*A[j][i];
        ck[i] = -A[j][i];

        fkx[i] = omega*alpha*L_phi_x[j][i];
        fky[i] = omega*alpha*L_phi_y[j][i];
      }

      tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fkx,ukx);
      tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fky,uky);

      for(int i=0;i<n-1;i++){
        fx[j][i] = ukx[i];
        fy[j][i] = uky[i];
      }
    }

    // Reapply periodicity boundary condition (for f) [confirmar se isto realmente é necessário]
    for(int j=1;j<m-1;j++){
      fx[j][n-1] = fx[j][0];
      fy[j][n-1] = fy[j][0];
    }

    // Solve for the deltas - step 2 (eta)
    for(int i=0;i<n-1;i++){
      // an[0] = -C[1][i];
      bn[0] = alpha + 2.*C[1][i];
      cn[0] = -C[1][i];

      fnx[0] = fx[1][i];// + C[1][i]*x[0][i];
      fny[0] = fy[1][i];// + C[1][i]*y[0][i];
      
      for(int j=2;j<m-2;j++){
        an[j-2] = -C[j][i];
        bn[j-1] = alpha + 2.*C[j][i];
        cn[j-1] = -C[j][i];

        fnx[j-1] = fx[j][i];
        fny[j-1] = fy[j][i];
        // printf("%i\n",j);
      }

      an[m-4] = -C[m-2][i];
      bn[m-3] = alpha + 2.*C[m-2][i];
      // cn[m-4] = -C[m-2][i];

      fnx[m-3] = fx[m-2][i];// + C[m-2][i]*x[m-1][i];
      fny[m-3] = fy[m-2][i];// + C[m-2][i]*y[m-1][i];

      tridiagonal_matrix_solver(m-2,an,bn,cn,fnx,unx);
      tridiagonal_matrix_solver(m-2,an,bn,cn,fny,uny);

      for(int j=1;j<m-1;j++){
        Delta_x[j][i] = unx[j-1];
        Delta_y[j][i] = uny[j-1];
      }
    }

    // Calculate new x and y
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        x[j][i] += Delta_x[j][i];
        y[j][i] += Delta_y[j][i];
      }
    }

    // Reapply periodicity boundary condition
    for(int j=1;j<m-1;j++){
      x[j][n-1] = x[j][0];
      y[j][n-1] = y[j][0];
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
  free(D);
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

/*
non-periodic version 
- gigiaero, 24/05/2026, 1739 hours
*/
void solve_adi_2d_rectangular_eom_np(int m,int n,double x[m][n],double y[m][n],
                                     sim_prmtrs *config,
                                     control_prmtrs *c_prmtrs){
  // Solver variables
  double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
  double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
  double (*Delta_x)[n] = calloc(m,sizeof *Delta_x);
  double (*Delta_y)[n] = calloc(m,sizeof *Delta_y);
  double (*A)[n] = calloc(m,sizeof *A);
  double (*B)[n] = calloc(m,sizeof *B);
  double (*C)[n] = calloc(m,sizeof *C);
  double (*D)[n] = calloc(m,sizeof *D);
  double (*fx)[n] = calloc(m,sizeof *fx);
  double (*fy)[n] = calloc(m,sizeof *fy);
  double *ak = malloc(sizeof(double)*(n-1));
  double *bk = malloc(sizeof(double)*(n-1));
  double *ck = malloc(sizeof(double)*(n-1));
  double *fkx = malloc(sizeof(double)*(n-1));
  double *fky = malloc(sizeof(double)*(n-1));
  double *ukx = malloc(sizeof(double)*(n-1));
  double *uky = malloc(sizeof(double)*(n-1));
  double *an = malloc(sizeof(double)*(m-3));
  double *bn = malloc(sizeof(double)*(m-2));
  double *cn = malloc(sizeof(double)*(m-3));
  double *fnx = malloc(sizeof(double)*(m-2));
  double *fny = malloc(sizeof(double)*(m-2));
  double *unx = malloc(sizeof(double)*(m-2));
  double *uny = malloc(sizeof(double)*(m-2));
  double omega = config->w;
  double alpha;
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
    alpha_sequence_aH(m,n,x,y,config,k);
    alpha_sequence(&alpha,&k,iter,config);

    calc_A(m,n,A,x,y);
    calc_B(m,n,B,x,y);
    calc_C(m,n,C,x,y);
    calc_D(m,n,D,x,y);

    // Calculate residual operators
    L_phi_eom(m,n,L_phi_x,L_phi_y,x,y,A,B,C,D,c_prmtrs);

    res_x = 0.;
    res_y = 0.;
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

    // Solve for the deltas - step 1 (ksi)
    for(int j=1;j<m-1;j++){
      bk[0] = alpha + 2.*A[j][1];
      ck[0] = -A[j][1];

      fkx[0] = omega*alpha*L_phi_x[j][1];
      fky[0] = omega*alpha*L_phi_y[j][1];

      for(int i=2;i<n-2;i++){
        ak[i-2] = -A[j][i];
        bk[i-1] = alpha + 2.*A[j][i];
        ck[i-1] = -A[j][i];

        fkx[i-1] = omega*alpha*L_phi_x[j][i];
        fky[i-1] = omega*alpha*L_phi_y[j][i];
      }

      ak[n-4] = -A[j][n-2];
      bk[n-3] = alpha + 2.*A[j][n-2];

      fkx[n-3] = omega*alpha*L_phi_x[j][n-2];
      fky[n-3] = omega*alpha*L_phi_y[j][n-2];

      tridiagonal_matrix_solver(n-2,ak,bk,ck,fkx,ukx);
      tridiagonal_matrix_solver(n-2,ak,bk,ck,fky,uky);

      for(int i=1;i<n-1;i++){
        fx[j][i] = ukx[i-1];
        fy[j][i] = uky[i-1];
      }
    }

    // // Reapply periodicity boundary condition (for f) [confirmar se isto realmente é necessário]
    // for(int j=1;j<m-1;j++){
    //   fx[j][n-1] = fx[j][0];
    //   fy[j][n-1] = fy[j][0];
    // }

    // Solve for the deltas - step 2 (eta)
    for(int i=1;i<n-1;i++){
      // an[0] = -C[1][i];
      bn[0] = alpha + 2.*C[1][i];
      cn[0] = -C[1][i];

      fnx[0] = fx[1][i];// + C[1][i]*x[0][i];
      fny[0] = fy[1][i];// + C[1][i]*y[0][i];
      
      for(int j=2;j<m-2;j++){
        an[j-2] = -C[j][i];
        bn[j-1] = alpha + 2.*C[j][i];
        cn[j-1] = -C[j][i];

        fnx[j-1] = fx[j][i];
        fny[j-1] = fy[j][i];
        // printf("%i\n",j);
      }

      an[m-4] = -C[m-2][i];
      bn[m-3] = alpha + 2.*C[m-2][i];
      // cn[m-4] = -C[m-2][i];

      fnx[m-3] = fx[m-2][i];// + C[m-2][i]*x[m-1][i];
      fny[m-3] = fy[m-2][i];// + C[m-2][i]*y[m-1][i];

      tridiagonal_matrix_solver(m-2,an,bn,cn,fnx,unx);
      tridiagonal_matrix_solver(m-2,an,bn,cn,fny,uny);

      for(int j=1;j<m-1;j++){
        Delta_x[j][i] = unx[j-1];
        Delta_y[j][i] = uny[j-1];
      }
    }

    // Calculate new x and y
    for(int j=1;j<m-1;j++){
      for(int i=1;i<n-1;i++){
        x[j][i] += Delta_x[j][i];
        y[j][i] += Delta_y[j][i];
      }
    }

    // // Reapply periodicity boundary condition
    // for(int j=1;j<m-1;j++){
    //   x[j][n-1] = x[j][0];
    //   y[j][n-1] = y[j][0];
    // }

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
  free(D);
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

/*
- gigiaero, 16/05/2026, 1223 hours
*/
void solve_af2_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config,control_prmtrs *c_prmtrs){
  // Solver variables
  double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
  double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
  double (*Delta_x)[n] = calloc(m,sizeof *Delta_x);
  double (*Delta_y)[n] = calloc(m,sizeof *Delta_y);
  double (*A)[n] = calloc(m,sizeof *A);
  double (*B)[n] = calloc(m,sizeof *B);
  double (*C)[n] = calloc(m,sizeof *C);
  double (*D)[n] = calloc(m,sizeof *D);
  double (*fx)[n] = calloc(m,sizeof *fx);
  double (*fy)[n] = calloc(m,sizeof *fy);
  double *ak = malloc(sizeof(double)*(n-1));
  double *bk = malloc(sizeof(double)*(n-1));
  double *ck = malloc(sizeof(double)*(n-1));
  double *fkx = malloc(sizeof(double)*(n-1));
  double *fky = malloc(sizeof(double)*(n-1));
  double *ukx = malloc(sizeof(double)*(n-1));
  double *uky = malloc(sizeof(double)*(n-1));
  double *an = malloc(sizeof(double)*(m-3));
  double *bn = malloc(sizeof(double)*(m-2));
  double *cn = malloc(sizeof(double)*(m-3));
  double *fnx = malloc(sizeof(double)*(m-2));
  double *fny = malloc(sizeof(double)*(m-2));
  double *unx = malloc(sizeof(double)*(m-2));
  double *uny = malloc(sizeof(double)*(m-2));
  double omega = config->w;
  double alpha;
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
    alpha_sequence_aH(m,n,x,y,config,k);
    alpha_sequence(&alpha,&k,iter,config);
    
    calc_A(m,n,A,x,y);
    calc_B(m,n,B,x,y);
    calc_C(m,n,C,x,y);
    calc_D(m,n,D,x,y);

    // Calculate residual operators
    L_phi_eom(m,n,L_phi_x,L_phi_y,x,y,A,B,C,D,c_prmtrs);

    res_x = 0.;
    res_y = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        if(fabs(L_phi_x[j][i]) > res_x)
          res_x = fabs(L_phi_x[j][i]);

        if(fabs(L_phi_y[j][i]) > res_y)
          res_y = fabs(L_phi_y[j][i]);
      }
    }

    printf("AF2 Iteration %010d | Res x %.6e | Res y %.6e\n",iter,res_x,res_y);

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
    for(int i=0;i<n-1;i++){
      bn[0] = alpha + C[1][i];
      cn[0] = -C[1][i];

      fnx[0] = alpha*omega*L_phi_x[1][i];
      fny[0] = alpha*omega*L_phi_y[1][i];
      
      for(int j=2;j<m-2;j++){
        an[j-2] = 0.;
        bn[j-1] = alpha + C[j][i];
        cn[j-1] = -C[j][i];

        fnx[j-1] = alpha*omega*L_phi_x[j][i];
        fny[j-1] = alpha*omega*L_phi_y[j][i];
      }

      an[m-4] = 0.;
      bn[m-3] = alpha + C[m-2][i];

      fnx[m-3] = alpha*omega*L_phi_x[m-2][i];// + C[m-2][i]*fx[m-1][i];
      fny[m-3] = alpha*omega*L_phi_y[m-2][i];

      tridiagonal_matrix_solver(m-2,an,bn,cn,fnx,unx);
      tridiagonal_matrix_solver(m-2,an,bn,cn,fny,uny);
      
      for(int j=1;j<m-1;j++){
        fx[j][i] = unx[j-1];
        fy[j][i] = uny[j-1];
      }
    }

    // Solve for the deltas - step 2 (ksi)
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        ak[i] = -A[j][i];
        bk[i] = 1. + 2.*A[j][i];
        ck[i] = -A[j][i];

        fkx[i] = fx[j][i] + Delta_x[j-1][i];
        fky[i] = fy[j][i] + Delta_y[j-1][i];
      }

      tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fkx,ukx);
      tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fky,uky);

      for(int i=0;i<n-1;i++){
        Delta_x[j][i] = ukx[i];
        Delta_y[j][i] = uky[i];
      }
    }

    // Calculate new x and y
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        x[j][i] += Delta_x[j][i];
        y[j][i] += Delta_y[j][i];
      }
    }

    // Reapply periodicity boundary condition
    for(int j=1;j<m-1;j++){
      x[j][n-1] = x[j][0];
      y[j][n-1] = y[j][0];
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
  free(D);
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

/*
line systems (which, apparently, actually works)
- gigiaero, 10/05/2026
*/
void solve_slor_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                   sim_prmtrs *config,control_prmtrs *c_prmtrs){
  // Solver variables
  double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
  double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
  double (*Delta_x)[n] = calloc(m,sizeof *Delta_x);
  double (*Delta_y)[n] = calloc(m,sizeof *Delta_y);
  double (*A)[n] = calloc(m,sizeof *A);
  double (*B)[n] = calloc(m,sizeof *B);
  double (*C)[n] = calloc(m,sizeof *C);
  double (*D)[n] = calloc(m,sizeof *D);
  double *ak = malloc(sizeof(double)*(n-1));
  double *bk = malloc(sizeof(double)*(n-1));
  double *ck = malloc(sizeof(double)*(n-1));
  double *fkx = malloc(sizeof(double)*(n-1));
  double *fky = malloc(sizeof(double)*(n-1));
  double *ukx = malloc(sizeof(double)*(n-1));
  double *uky = malloc(sizeof(double)*(n-1));
  double omega = config->w;
  double inv_r = 1./config->r;
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
    calc_A(m,n,A,x,y);
    calc_B(m,n,B,x,y);
    calc_C(m,n,C,x,y);
    calc_D(m,n,D,x,y);

    // Calculate residual operators
    L_phi_eom(m,n,L_phi_x,L_phi_y,x,y,A,B,C,D,c_prmtrs);

    res_x = 0.;
    res_y = 0.;
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        if(fabs(L_phi_x[j][i]) > res_x)
          res_x = fabs(L_phi_x[j][i]);

        if(fabs(L_phi_y[j][i]) > res_y)
          res_y = fabs(L_phi_y[j][i]);
      }
    }

    printf("SLOR Iteration %010d | Res x %.6e | Res y %.6e\n",iter,res_x,res_y);

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

    // Calculate deltas
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        ak[i] = A[j][i];
        bk[i] = -2.*(inv_r + A[j][i]);
        ck[i] = A[j][i];

        fkx[i] = -omega*L_phi_x[j][i] - Delta_x[j-1][i];
        fky[i] = -omega*L_phi_y[j][i] - Delta_y[j-1][i];
      }

      tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fkx,ukx);
      tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fky,uky);

      for(int i=0;i<n-1;i++){
        Delta_x[j][i] = ukx[i];
        Delta_y[j][i] = uky[i];
      }
    }

    // Calculate new x and y
    for(int j=1;j<m-1;j++){
      for(int i=0;i<n-1;i++){
        x[j][i] += Delta_x[j][i];
        y[j][i] += Delta_y[j][i];
      }
    }

    // Reapply periodicity boundary condition
    for(int j=1;j<m-1;j++){
      x[j][n-1] = x[j][0];
      y[j][n-1] = y[j][0];
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
  free(D);
  free(ak);
  free(bk);
  free(ck);
  free(ukx);
  free(uky);
  free(fkx);
  free(fky);
  free(filename_save_x);
  free(filename_save_y);
  free(buffer);
  free(filename_log_x);
  free(filename_log_y);
  fclose(file_log_x);
  fclose(file_log_y);
}

// /*
// original version (column systems)
// - gigiaero, 10/05/2026, 1928 hours
// */
// void solve_slor_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
//                                    sim_prmtrs *config,control_prmtrs *c_prmtrs){
//   // Solver variables
//   double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
//   double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
//   double (*Delta_x)[n] = calloc(m,sizeof *Delta_x);
//   double (*Delta_y)[n] = calloc(m,sizeof *Delta_y);
//   double (*A)[n] = calloc(m,sizeof *A);
//   double (*B)[n] = calloc(m,sizeof *B);
//   double (*C)[n] = calloc(m,sizeof *C);
//   double (*D)[n] = calloc(m,sizeof *D);
//   // double (*fx)[n] = calloc(m,sizeof *fx);
//   // double (*fy)[n] = calloc(m,sizeof *fy);
//   // double *ak = malloc(sizeof(double)*(n-1));
//   // double *bk = malloc(sizeof(double)*(n-1));
//   // double *ck = malloc(sizeof(double)*(n-1));
//   // double *fkx = malloc(sizeof(double)*(n-1));
//   // double *fky = malloc(sizeof(double)*(n-1));
//   // double *ukx = malloc(sizeof(double)*(n-1));
//   // double *uky = malloc(sizeof(double)*(n-1));
//   double *an = malloc(sizeof(double)*(m-3));
//   double *bn = malloc(sizeof(double)*(m-2));
//   double *cn = malloc(sizeof(double)*(m-3));
//   double *fnx = malloc(sizeof(double)*(m-2));
//   double *fny = malloc(sizeof(double)*(m-2));
//   double *unx = malloc(sizeof(double)*(m-2));
//   double *uny = malloc(sizeof(double)*(m-2));
//   double omega = config->w;
//   double inv_r = 1./config->r;
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
//     calc_A(m,n,A,x,y);
//     calc_B(m,n,B,x,y);
//     calc_C(m,n,C,x,y);
//     calc_D(m,n,D,x,y);

//     // Calculate residual operators
//     L_phi_eom(m,n,L_phi_x,L_phi_y,x,y,A,B,C,D,c_prmtrs);

//     res_x = 0.;
//     res_y = 0.;
//     for(int j=1;j<m-1;j++){
//       for(int i=0;i<n-1;i++){
//         if(fabs(L_phi_x[j][i]) > res_x)
//           res_x = fabs(L_phi_x[j][i]);

//         if(fabs(L_phi_y[j][i]) > res_y)
//           res_y = fabs(L_phi_y[j][i]);
//       }
//     }

//     printf("SLOR Iteration %010d | Res x %.6e | Res y %.6e\n",iter,res_x,res_y);

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

//     // Solve for the deltas
//     // Column i = 0
//     bn[0] = -2.*(inv_r + C[1][0]);
//     cn[0] = C[1][0];

//     fnx[0] = -omega*L_phi_x[1][0] - Delta_x[1][n-2] - C[1][0]*Delta_x[0][0]*0;
//     fny[0] = -omega*L_phi_y[1][0] - Delta_y[1][n-2] - C[1][0]*Delta_y[0][0]*0;
    
//     for(int j=2;j<m-2;j++){
//       an[j-2] = C[j][0];
//       bn[j-1] = -2.*(inv_r + C[j][0]);
//       cn[j-1] = C[j][0];

//       fnx[j-1] = -omega*L_phi_x[j][0] - Delta_x[j][n-2];
//       fny[j-1] = -omega*L_phi_y[j][0] - Delta_y[j][n-2];
//     }

//     an[m-4] = C[m-2][0];
//     bn[m-3] = -2.*(inv_r + C[m-2][0]);

//     fnx[m-3] = -omega*L_phi_x[m-2][0] - Delta_x[m-2][n-2] - C[m-2][0]*Delta_x[m-1][0]*0;
//     fny[m-3] = -omega*L_phi_y[m-2][0] - Delta_y[m-2][n-2] - C[m-2][0]*Delta_y[m-1][0]*0;

//     tridiagonal_matrix_solver(m-2,an,bn,cn,fnx,unx);
//     tridiagonal_matrix_solver(m-2,an,bn,cn,fny,uny);

//     for(int j=1;j<m-1;j++){
//       Delta_x[j][0] = unx[j-1];
//       Delta_y[j][0] = uny[j-1];
//     }

//     // Rest of the domain
//     for(int i=1;i<n-1;i++){
//       bn[0] = -2.*(inv_r + C[1][i]);
//       cn[0] = C[1][i];

//       fnx[0] = -omega*L_phi_x[1][i] - Delta_x[1][i-1] - C[1][i]*Delta_x[0][i]*0;
//       fny[0] = -omega*L_phi_y[1][i] - Delta_y[1][i-1] - C[1][i]*Delta_y[0][i]*0;
      
//       for(int j=2;j<m-2;j++){
//         an[j-2] = C[j][i];
//         bn[j-1] = -2.*(inv_r + C[j][i]);
//         cn[j-1] = C[j][i];

//         fnx[j-1] = -omega*L_phi_x[j][i] - Delta_x[j][i-1];
//         fny[j-1] = -omega*L_phi_y[j][i] - Delta_y[j][i-1];
//       }

//       an[m-4] = C[m-2][i];
//       bn[m-3] = -2.*(inv_r + C[m-2][i]);

//       fnx[m-3] = -omega*L_phi_x[m-2][i] - Delta_x[m-2][i-1] - C[m-2][i]*Delta_x[m-1][i]*0;
//       fny[m-3] = -omega*L_phi_y[m-2][i] - Delta_y[m-2][i-1] - C[m-2][i]*Delta_y[m-1][i]*0;

//       tridiagonal_matrix_solver(m-2,an,bn,cn,fnx,unx);
//       tridiagonal_matrix_solver(m-2,an,bn,cn,fny,uny);

//       for(int j=1;j<m-1;j++){
//         Delta_x[j][i] = unx[j-1];
//         Delta_y[j][i] = uny[j-1];
//       }
//     }

//     // Calculate new x and y
//     for(int j=1;j<m-1;j++){
//       for(int i=0;i<n-1;i++){
//         x[j][i] += Delta_x[j][i];
//         y[j][i] += Delta_y[j][i];
//       }
//     }

//     // Reapply periodicity boundary condition
//     for(int j=1;j<m-1;j++){
//       x[j][n-1] = x[j][0];
//       y[j][n-1] = y[j][0];
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
//   free(D);
//   // free(fx);
//   // free(fy);
//   // free(ak);
//   // free(bk);
//   // free(ck);
//   // free(ukx);
//   // free(uky);
//   // free(fkx);
//   // free(fky);
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
//   fclose(file_log_x);
// }

double uniform_scheme_der1_o2_central_prdc_ksi(int m,int n,double phi[m][n],
                                               int j){
  return (phi[j][1] - phi[j][n-2])*.5;
}

double uniform_scheme_der2_o2_central_prdc_ksi(int m,int n,double phi[m][n],
                                               int j,int axis){
  switch(axis){
    case 1: // Horizontal
      return phi[j][1] - 2.*phi[j][0] + phi[j][n-2];
      break;

    case 3: // Mixed
      return (phi[j+1][1] - phi[j-1][1] - phi[j+1][n-2] + phi[j-1][n-2])*.25;
      break;

    default:
      puts("uniform_scheme_der2_o2_central_prdc_ksi: invalid axis");
      exit(15);
  }
}