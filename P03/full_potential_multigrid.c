#include<stdio.h>
#include<stdlib.h>
#include"../include/fullp_lib.h"

// A FAZER: funções que salvam em arquivos de texto os parâmetros da malha e da 
// solução

int main(){
  /* MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* Biconvex airfoil, 10% thickness */
  char *fname_msh_x[] = {"../P02/eom_m81_n241_biair_t01/eom_x_iter_0000000065.dat",
                         "../P02/eom_m41_n121_biair_t01/eom_x_iter_0000000064.dat",
                         "../P02/eom_m21_n61_biair_t01/eom_x_iter_0000000047.dat",
                         "../P02/eom_m11_n31_biair_t01/eom_x_iter_0000000168.dat"};
  char *fname_msh_y[] = {"../P02/eom_m81_n241_biair_t01/eom_y_iter_0000000065.dat",
                         "../P02/eom_m41_n121_biair_t01/eom_y_iter_0000000064.dat",
                         "../P02/eom_m21_n61_biair_t01/eom_y_iter_0000000047.dat",
                         "../P02/eom_m11_n31_biair_t01/eom_y_iter_0000000168.dat"};
  int n[] = {241,121,61,31};
  int m[] = {81,41,21,11};

  /* FLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  sim_prmtrs config;
  config.Ntype = 3;
  config.w = 2.;
  config.r = 1.;
  config.alpha_seq = 0;
  config.alpha = 1;
  config.alpha_H = 1e3;
  config.set_alpha_H = 0;
  config.M = 5;
  config.max_iter = 1000;
  config.qtimes = 1;
  config.save_i_c = 1;
  config.save_last_only = 1;
  config.eps = 1.e-6;
  // config.eps = 1.e-8;
  char output_file[] = "results/fullp";

  fullp_prmtrs fp_prmtrs;
  fp_prmtrs.alpha = 0.;
  fp_prmtrs.Ma = .7;
  fp_prmtrs.C = 1.;
  fp_prmtrs.beta_sub = .3;
  fp_prmtrs.beta_super = 4.5;
  fp_prmtrs.lift = 0;
  fp_prmtrs.Rg = .1;

  /* FLOW (solution) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  config.casename = malloc(sizeof(char)*200);
  sprintf(config.casename,"%s",output_file);
  fullp_multigrid(m,n,&config,&fp_prmtrs,fname_msh_x,fname_msh_y);

  free(config.casename);

  return 0;
}