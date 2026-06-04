#include<stdio.h>
#include<stdlib.h>
#include"../include/eom_lib.h"
#include"../include/fullp_lib.h"
#include"../include/mesh.h"

// nome padrão dos arquivos de malha: fullp_mesh_x.dat e fullp_mesh_y.dat
// init_only = 1 -> renomear a malha inicial pro nome padrão
// init_only = 0 -> renomear a malha final pro nome padrão

// A FAZER: funções que salvam em arquivos de texto os parâmetros da malha e da 
// solução

int main(){
  /* MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  char filename_msh_x[] = "../P02/results/eom_x_iter_0000000047.dat";
  char filename_msh_y[] = "../P02/results/eom_y_iter_0000000047.dat";
  int m = 15;
  int n = 93;

  /* FLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  sim_prmtrs config;
  config.Ntype = 2;
  config.w = 2.;
  config.r = 1.;
  config.alpha_seq = 0;
  config.alpha = 1;
  config.alpha_H = 1e-1;
  config.set_alpha_H = 0;
  config.M = 5;
  config.max_iter = 500000;
  config.qtimes = 72605;
  config.save_i_c = 1;
  config.save_last_only = 1;
  config.eps = 1.e-6;
  char output_file[] = "results/fullp";

  fullp_prmtrs fp_prmtrs;
  fp_prmtrs.alpha = 0.;
  fp_prmtrs.Ma = .7;
  fp_prmtrs.C = 1.;
  fp_prmtrs.lift = 0;

  /* FLOW (solution) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  config.casename = malloc(sizeof(char)*200);
  sprintf(config.casename,"%s",output_file);
  evaluate_delta_form_fullp(m,n,&config,&fp_prmtrs,filename_msh_x,
                            filename_msh_y);

  free(config.casename);

  return 0;
}