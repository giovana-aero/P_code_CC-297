#include<math.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/eom_lib.h"
#include"../include/num_methods.h"

void calc_xy_eta(double *xy_eta,int m,int n,double x[m][n],double y[m][n],
                 int i,int j,double S_eta){
  double x_ksi;
  double y_ksi;

  if(i == 0){
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

void parabolic_mesh(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh){
  double x_int,y_int,x0,y0;
  double delta_s;
  double R;
  double *eps_switch = malloc(sizeof(double)*m);
  double *s = malloc(sizeof(double)*m);
  double *xy_eta = malloc(sizeof(double)*2);
  double S_eta;

  linspace(eps_switch,0.,1.,m);
  // cosspace(s,0.,1.,m,1);
  linspace(s,0.,1.,m);

  // for(int j=1;j<m-3;j++){
    for(int j=1;j<m-1;j++){
    for(int i=0;i<n-1;i++){
      // Local reference grid at j
      calc_xy_eta(xy_eta,m,n,x,y,i,j,S_eta);
      // x0 = x[j-1][i] + xy_eta[0];
      // y0 = y[j-1][i] + xy_eta[1];
      x0 = x[j-1][i] + uniform_scheme_der1_o1_forward(m,n,x,i,j,2);
      y0 = y[j-1][i] + uniform_scheme_der1_o1_forward(m,n,y,i,j,2);

      R = ((double) m - j);
      delta_s = (s[j] - s[j-1])/(s[m-1] - s[j-1]);
      S_eta = delta_s*R;

      x_int = x[j-1][i] + delta_s*(x[m-1][i] - x[j-1][i]);
      y_int = y[j-1][i] + delta_s*(y[m-1][i] - y[j-1][i]);

      x[j][i] = eps_switch[j]*x_int + (1. - eps_switch[j])*x0;
      y[j][i] = eps_switch[j]*y_int + (1. - eps_switch[j])*y0; 

      // Local reference grid at j+1
    }

    x[j][n-1] = x[j][0];
    y[j][n-1] = y[j][0];

    // break;
    // Solve for actual values at j

  }

  /* adicionar aqui a solução pra penúltima linha */

  free(eps_switch);
  free(xy_eta);
  free(s);
}
