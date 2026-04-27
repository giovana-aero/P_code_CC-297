#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"../include/bi_air_lib.h"
#include"../include/eom_lib.h"

#define pi 3.1415926535897932384626433

/*
- gigiaero, 27/04/2026, 1310 hours
*/
void calc_A(int m,int n,double A[m][n],double x[m][n],double y[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      A[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,2),2) + 
                pow(uniform_scheme_der1_o2_central(m,n,y,i,j,2),2);
  }
}

/*
- gigiaero, 27/04/2026, 1310 hours
*/
void calc_B(int m,int n,double B[m][n],double x[m][n],double y[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      B[j][i] = uniform_scheme_der1_o2_central(m,n,x,i,j,1)*
                uniform_scheme_der1_o2_central(m,n,x,i,j,2) + 
                uniform_scheme_der1_o2_central(m,n,y,i,j,1)*
                uniform_scheme_der1_o2_central(m,n,y,i,j,2);
  }
}

/*
- gigiaero, 27/04/2026, 1310 hours
*/
void calc_C(int m,int n,double C[m][n],double x[m][n],double y[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      C[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,1),2) + 
                pow(uniform_scheme_der1_o2_central(m,n,y,i,j,1),2);
  }
}

/*
- gigiaero, 27/04/2026, 1310 hours
*/
void calc_D(int m,int n,double D[m][n],double x[m][n],double y[m][n]){
  for(int j=0;j<m;j++){
    for(int i=0;i<n;i++)
      D[j][i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,1)*
                    uniform_scheme_der1_o2_central(m,n,y,i,j,2) - 
                    uniform_scheme_der1_o2_central(m,n,x,i,j,2)*
                    uniform_scheme_der1_o2_central(m,n,y,i,j,1),2);
  }
}

/*
- gigiaero, 24/04/2026, 1352 hours
*/
void cosspace(double *x,double xi,double xf,int n,int half){
  double midp;
  double th_i = pi/(((double) n) - 1.);
  double th = th_i;

  if(half){
    midp = xf - xi;
    th_i /= 2.;
  }
  else 
    midp = (xf - xi)/2.;

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
    // prmtrs[i+1] = v_ex[i];
    v_in[i] = prmtrs[i+1+n+2];
    // prmtrs[i+1+n+2] = v_in[i];
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
      K = factorial(n)/(factorial(r)*factorial(n-r));
      sum1 += A1[r]*K*pow(x[i],r)*pow(1. - x[i],n - r);
      sum2 += A2[r]*K*pow(x[i],r)*pow(1. - x[i],n - r);
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
*/
void cst_prmtrs(int af,double *prmtrs){
  int n = prmtrs[0];
  double *v_ex = malloc(sizeof(double)*(n+2));
  double *v_in = malloc(sizeof(double)*(n+2));
  char *filename = malloc(sizeof(char)*100);
  sprintf(filename,"./reverse_cst/");
  FILE *input;
  
  switch(af){
    case 1: // whitcomb supercritical, n = 10
      strcat(filename,"whitcomb_cst_prmtrs.dat");
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
    prmtrs[i+1] = v_ex[i];
    prmtrs[i+1+n+2] = v_in[i];
  }

  // print_1d_array((n+2)*2 + 1,prmtrs);

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
    linspace(th,2*pi,0.,n);
  else
    linspace(th,0.,2*pi,n);

  for(int i=0;i<n;i++){
    x[i] = prmtrs[2] + prmtrs[0]*cos(th[i]);
    y[i] = prmtrs[3] + prmtrs[1]*sin(th[i]);
  }
    
  free(th);
}

/*
- gigiaero, 26/04/2026, 1003 hours
*/
void evaluate_delta_form_eom(sim_prmtrs *config,msh_prmtrs *msh){
  double (*x)[msh->IMAX] = calloc(msh->JMAX,sizeof *x);
  double (*y)[msh->IMAX] = calloc(msh->JMAX,sizeof *y);

  initialize_mesh(msh->JMAX,msh->IMAX,x,y,msh);

  switch(config->Ntype){
    case 1:
      // slor
      break;
    
    case 2:
      // adi
      break;
      
    default:
      puts("evaluate_delta_form_eom: Invalid Ntype");
      exit(32);
  }

  

  // print_2d_array_to_file(msh.JMAX,msh.IMAX,x,"mesh_x.dat",0);
  // print_2d_array_to_file(msh.JMAX,msh.IMAX,y,"mesh_y.dat",0);


  free(x);
  free(y);
}

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
*/
void init_type2(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  double *vx = malloc(sizeof(double)*msh->JMAX);
  double *vy = malloc(sizeof(double)*msh->JMAX);

  for(int i=0;i<msh->IMAX;i++){
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
*/
void init_type3(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  double *vx = malloc(sizeof(double)*msh->JMAX);
  double *vy = malloc(sizeof(double)*msh->JMAX);

  for(int i=0;i<msh->IMAX;i++){
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
      puts("init_type 4 pending!");
      exit(4);
      break;

    default:
      puts("initialize_mesh: invalid init_type");
      exit(4);
  }

  free(x_axis);
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

void solve_adi_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config){
  // Solver variables
  double (*L_phi_x)[n] = calloc(m,sizeof *L_phi_x);
  double (*L_phi_y)[n] = calloc(m,sizeof *L_phi_y);
  double (*A)[n] = calloc(m,sizeof *A);
  double (*B)[n] = calloc(m,sizeof *B);
  double (*C)[n] = calloc(m,sizeof *C);
  double (*D)[n] = calloc(m,sizeof *D);
  int iter = 0;
  // Save files
  char *filename_save_x = malloc(sizeof(char)*200);
  char *filename_save_y = malloc(sizeof(char)*200);
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
  sprintf(filename_save_x,"%s_x_iter_",config->casename);
  sprintf(filename_save_y,"%s_y_iter_",config->casename);
  find_str_end(filename_save_x,&str_end_idx);

  for(iter;iter<=config->max_iter;iter++){
    // calcular operadores residuo

    // teste de convergencia
    // break;

    // calculo das correcoes

    // atualizacao das coordenadas


  }

  free(L_phi_x);
  free(L_phi_y);
  free(A);
  free(B);
  free(C);
  free(D);
  free(filename_save_x);
  free(filename_save_y);
  free(buffer);
  free(filename_log);
  fclose(file_log);
}

/*
periodic tridiagonal matrices
- gigiaero, 26/04/2026, 2236 hours
*/
void tridiagonal_pmatrix_solver(int n,double *a,double *b,double *c,double *f,
                                double *u){
  double *x = malloc(sizeof(double)*(n-1));
  double *w = malloc(sizeof(double)*n);
  double *k = malloc(sizeof(double)*(n-1));
  double *g = malloc(sizeof(double)*(n-2));
  double *h = malloc(sizeof(double)*(n-2));
  double *y = malloc(sizeof(double)*n);

  // Obtain matrix coefficients
  h[0] = c[n-1];
  w[0] = b[0];
  g[0] = a[0]/w[0];
  x[0] = c[0]/w[0];
  k[0] = a[1];

  for(int i=1;i<n-2;i++){
    h[i] = -h[i-1]*x[i-1];
    w[i] = b[i] - k[i-1]*x[i-1];
    g[i] = -g[i-1]*k[i-1]/w[i];
    x[i] = c[i]/w[i];
    k[i] = a[i+1];
  }

  w[n-2] = b[n-2] - k[n-3]*x[n-3];
  x[n-2] = (c[n-2] - g[n-3]*k[n-3])/w[n-2];
  k[n-2] = a[n-1] - h[n-3]*x[n-3];

  w[n-1] = b[n-1] - k[n-2]*x[n-2];
  for(int i=0;i<n-2;i++)
    w[n-1] -= g[i]*h[i];
  
  // Calculate y
  y[0] = f[0]/w[0];
  y[n-1] = f[n-1] - h[0]*y[0];

  for(int i=1;i<n-1;i++){
    y[i] = (f[i] - k[i-1]*y[i-1])/w[i];
    
    if(i != n-2)
      y[n-1] -= h[i]*y[i];

    else
      y[n-1] -= k[i]*y[i];
  }

  y[n-1] /= w[n-1];
  
  // Calculate solution
  u[n-1] = y[n-1];
  u[n-2] = y[n-2] - u[n-1]*x[n-2];
  for(int i=n-3;i>=0;i--)
    u[i] = y[i] - u[i+1]*x[i] - u[n-1]*g[i];

  free(x);
  free(w);
  free(k);
  free(g);
  free(h);
  free(y);
}