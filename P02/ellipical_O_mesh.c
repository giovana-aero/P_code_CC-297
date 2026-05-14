#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/bi_air_lib.h"
#include"../include/eom_lib.h"
#include"../include/mesh.h"

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

int main(){
  // Solution configurations
  sim_prmtrs config;
  config.Ntype = 1;
  config.w = 2.;
  // config.w = 1.;
  config.r = 1.6;
  // config.r = 1;
  config.alpha = .01;
  // config.alpha = 1.;
  config.max_iter = 100000;
  config.qtimes = 25000;
  config.save_i_c = 1;
  config.save_last_only = 1;
  config.eps = 1.e-6; // Convergence criterion
  char output_file[] = "results/eom";

  // Mesh parameters
  msh_prmtrs msh;
  msh.end_prmtrs = malloc(sizeof(double)*4);
  /* IMAX */
  msh.IMAX = 93;
  /* JMAX */
  msh.JMAX = 15;
  /* c */
  msh.c = 1.;
  /* end_prmtrs */
  msh.end_prmtrs[0] = 6.5*msh.c;
  msh.end_prmtrs[1] = 6.5*msh.c;
  msh.end_prmtrs[2] = msh.c*.5;
  msh.end_prmtrs[3] = 0.;
  /* init_type */
  msh.init_type = 4;
  int init_only = 1; // Initialize only, do not solve

  /* af_type */
  int n = 10; // cst - bernstein polynomial order
  msh.af_prmtrs = malloc(sizeof(double)*((n+2)*2 + 1));
  msh.af_type = 3;
  /* af_prmtrs (bi_air) */
  // msh.af_prmtrs[0] = 0.1;
  /* af_prmtrs (naca4) */
  // msh.af_prmtrs[0] = 8.;
  // msh.af_prmtrs[1] = 4.;
  // msh.af_prmtrs[2] = 12.;
  /* af_prmtrs (cst) */
  msh.af_prmtrs[0] = n;
  cst_prmtrs(1,msh.af_prmtrs);

  // P & Q control functions
  control_prmtrs c_prmtrs;
  c_prmtrs.L = 1;
  c_prmtrs.M = 1;
  malloc_c_prmtrs(&c_prmtrs);
  c_prmtrs.al[0] = 0;
  c_prmtrs.bm[0] = 0;
  c_prmtrs.cl[0] = 10;
  c_prmtrs.dm[0] = 1;
  c_prmtrs.ksi_l[0] = 20;
  c_prmtrs.eta_l[0] = 0;
  c_prmtrs.ksi_m[0] = 0;
  c_prmtrs.eta_m[0] = 0;

  // Solve
  config.casename = malloc(sizeof(char)*200);
  sprintf(config.casename,"%s",output_file);
  evaluate_delta_form_eom(&config,&msh,&c_prmtrs,init_only);

  free(msh.af_prmtrs);
  free(msh.end_prmtrs);
  free(config.casename);
  free(c_prmtrs.al);
  free(c_prmtrs.bm);
  free(c_prmtrs.cl);
  free(c_prmtrs.dm);
  free(c_prmtrs.ksi_l);
  free(c_prmtrs.ksi_m);
  free(c_prmtrs.eta_l);
  free(c_prmtrs.eta_m);

  // free(x_tmp);
  // free(y_tmp);

  return 0;
}