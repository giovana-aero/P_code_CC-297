#include"../include/bi_air_lib.h"
#include"../include/b_conditions.h"

/*
- gigiaero, 22/03/2026, 2131 hours
*/
void biconvex_airfoil_mesh(bi_air_mesh *b_a_m,double *x, double*y){

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