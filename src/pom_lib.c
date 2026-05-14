#include<math.h>
#include<stdlib.h>
// #include"../include/2d_arrays.h"
#include"../include/eom_lib.h"
#include"../include/num_methods.h"

/*
- gigiaero, 14/05/2026, 1646 hours
*/
void calc_A_line(int m,int n,double *A,double x[m][n],double y[m][n],int j){
  for(int i=0;i<n-1;i++)
    A[i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,2),2.) + 
           pow(uniform_scheme_der1_o2_central(m,n,y,i,j,2),2.);
  
}

/*
- gigiaero, 14/05/2026, 1646 hours
*/
void calc_B_line(int m,int n,double *B,double x[m][n],double y[m][n],int j){
  B[0] = uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j)*
         uniform_scheme_der1_o2_central(m,n,x,0,j,2) + 
         uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j)*
         uniform_scheme_der1_o2_central(m,n,y,0,j,2);

  for(int i=1;i<n-1;i++)
    B[i] = uniform_scheme_der1_o2_central(m,n,x,i,j,1)*
           uniform_scheme_der1_o2_central(m,n,x,i,j,2) + 
           uniform_scheme_der1_o2_central(m,n,y,i,j,1)*
           uniform_scheme_der1_o2_central(m,n,y,i,j,2);
  
}

/*
- gigiaero, 14/05/2026, 1646 hours
*/
void calc_C_line(int m,int n,double *C,double x[m][n],double y[m][n],int j){
  C[0] = pow(uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j),2.) + 
         pow(uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j),2.);

  for(int i=1;i<n-1;i++)
    C[i] = pow(uniform_scheme_der1_o2_central(m,n,x,i,j,1),2.) + 
           pow(uniform_scheme_der1_o2_central(m,n,y,i,j,1),2.);
  
}

/*
- gigiaero, 14/05/2026, 1733 hours 
*/
void calc_xy_eta(double *xy_eta,int m,int n,double x[m][n],double y[m][n],
                 int i,int j,double S_eta){
  double x_ksi;
  double y_ksi;

  if(i==0){
    x_ksi = uniform_scheme_der1_o2_central_prdc_ksi(m,n,x,j-1);
    y_ksi = uniform_scheme_der1_o2_central_prdc_ksi(m,n,y,j-1);
  }
  else{
    x_ksi = uniform_scheme_der1_o2_central(m,n,x,i,j-1,1);
    y_ksi = uniform_scheme_der1_o2_central(m,n,y,i,j-1,1);
  }

  xy_eta[0] = -y_ksi*S_eta/sqrt(x_ksi*x_ksi + y_ksi*y_ksi);
  xy_eta[1] = x_ksi*S_eta/sqrt(x_ksi*x_ksi + y_ksi*y_ksi);
}

/*
- gigiaero, 14/05/2026, 1733 hours 
*/
void exspace(double *x,double a,int n){
  for(int i=0;i<n;i++)
    x[i] = exp(((double) i)*a);
}

/*
- gigiaero, 14/05/2026, 1733 hours 
*/
void loc_ref_grid(int m,int n,double x[m][n],double y[m][n],double *eps_switch,
                  double *s,int j){
  double x_int,y_int,x0,y0;
  double S_eta,delta_s;
  double R;
  double *xy_eta = malloc(sizeof(double)*2);

  for(int i=0;i<n-1;i++){
    R = ((double) m - j);
    delta_s = (s[j] - s[j-1])/(s[m-1] - s[j-1]);
    S_eta = delta_s*R;

    calc_xy_eta(xy_eta,m,n,x,y,i,j,S_eta);
    x0 = x[j-1][i] + xy_eta[0];
    y0 = y[j-1][i] + xy_eta[1];
    // x0 = x[j-1][i] + uniform_scheme_der1_o1_forward(m,n,x,i,j,2);
    // y0 = y[j-1][i] + uniform_scheme_der1_o1_forward(m,n,y,i,j,2);

    x_int = x[j-1][i] + delta_s*(x[m-1][i] - x[j-1][i]);
    y_int = y[j-1][i] + delta_s*(y[m-1][i] - y[j-1][i]);

    x[j][i] = eps_switch[j]*x_int + (1. - eps_switch[j])*x0;
    y[j][i] = eps_switch[j]*y_int + (1. - eps_switch[j])*y0; 
  }

  x[j][n-1] = x[j][0];
  y[j][n-1] = y[j][0];

  free(xy_eta);
}

/*
- gigiaero, 14/05/2026, 1733 hours 
*/
void parabolic_mesh(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  double *eps_switch = malloc(sizeof(double)*m);
  double *s = malloc(sizeof(double)*m);
  double *A = malloc(sizeof(double)*n);
  double *B = malloc(sizeof(double)*n);
  double *C = malloc(sizeof(double)*n);
  double *ak = malloc(sizeof(double)*(n-1));
  double *bk = malloc(sizeof(double)*(n-1));
  double *ck = malloc(sizeof(double)*(n-1));
  double *fkx = malloc(sizeof(double)*(n-1));
  double *fky = malloc(sizeof(double)*(n-1));
  double *ukx = malloc(sizeof(double)*(n-1));
  double *uky = malloc(sizeof(double)*(n-1));

  linspace(eps_switch,0.,1.,m);
  // cosspace(s,0.,1.,m,1);
  // linspace(s,0.,1.,m);
  exspace(s,.15,m); // the only one that properly works every time, it seems.
                    // the .15 value was chosen kind of empirically.

  for(int j=1;j<m-1;j++){
    loc_ref_grid(m,n,x,y,eps_switch,s,j);
    
    if(j<m-2)
      loc_ref_grid(m,n,x,y,eps_switch,s,j+1);

    calc_A_line(m,n,A,x,y,j);
    calc_B_line(m,n,B,x,y,j);
    calc_C_line(m,n,C,x,y,j);
    
    for(int i=0;i<n-1;i++){
      ak[i] = 2.*A[i];
      bk[i] = -4.*(A[i] + C[i]);
      ck[i] = 2.*A[i];

      fkx[i] = B[i]*(x[j+1][i+1] - x[j+1][i-1] - x[j-1][i+1] + x[j-1][i-1])
               - 2.*C[i]*(x[j+1][i] + x[j-1][i]);
      fky[i] = B[i]*(y[j+1][i+1] - y[j+1][i-1] - y[j-1][i+1] + y[j-1][i-1])
               - 2.*C[i]*(y[j+1][i] + y[j-1][i]);
    }

    tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fkx,ukx);
    tridiagonal_pmatrix_solver(n-1,ak,bk,ck,fky,uky);

    for(int i=0;i<n-1;i++){
      x[j][i] = ukx[i];
      y[j][i] = uky[i];
    }

    x[j][n-1] = x[j][0];
    y[j][n-1] = y[j][0];
  }

  free(eps_switch);
  free(s);
  free(A);
  free(B);
  free(C);
  free(ak);
  free(bk);
  free(ck);
  free(fkx);
  free(fky);
  free(ukx);
  free(uky);
}
