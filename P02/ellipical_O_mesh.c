#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/bi_air_lib.h"
#include"../include/eom_lib.h"
#include"../include/mesh.h"

int main(){
  // Solution configurations
  sim_prmtrs config;
  config.Ntype = 5;
  config.r = 1.6; //slor
  // config.w = 1.5;
  config.max_iter = 40;
  config.qtimes = 1;
  config.save_i_c = 1;
  config.save_last_only = 0;
  config.eps = 1.e-6; // Convergence criterion
  char output_file[] = "results/eom";

  msh_prmtrs msh;
  int n = 10; // cst - bernstein polynomial order
  msh.af_prmtrs = malloc(sizeof(double)*((n+2)*2 + 1));
  msh.end_prmtrs = malloc(sizeof(double)*4);
  /* IMAX */
  msh.IMAX = 93;
  /* JMAX */
  msh.JMAX = 15;
  /* c */
  msh.c = 1.;
  /* af_type */
  msh.af_type = 1;
  /* af_prmtrs (bi_air) */
  msh.af_prmtrs[0] = 0.1;
  /* af_prmtrs (naca4) */
  // msh.af_prmtrs[0] = 8.;
  // msh.af_prmtrs[1] = 4.;
  // msh.af_prmtrs[2] = 12;
  /* af_prmtrs (cst) */
  // msh.af_prmtrs[0] = n;
  cst_prmtrs(1,msh.af_prmtrs);
  /* end_prmtrs */
  msh.end_prmtrs[0] = 6.5*msh.c;
  msh.end_prmtrs[1] = 6.5*msh.c;
  msh.end_prmtrs[2] = msh.c*.5;
  msh.end_prmtrs[3] = 0.;
  /* init_type */
  msh.init_type = 1;



  

  


  


  free(msh.af_prmtrs);
  free(msh.end_prmtrs);


  // free(x_tmp);
  // free(y_tmp);

  return 0;
}