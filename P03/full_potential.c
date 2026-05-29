#include<stdio.h>
#include<stdlib.h>
#include"../include/eom_lib.h"
#include"../include/mesh.h"

// nome padrão dos arquivos de malha: fullp_mesh_x.dat e fullp_mesh_y.dat
// init_only = 1 -> renomear a malha inicial pro nome padrão
// init_only = 0 -> renomear a malha final pro nome padrão

// A FAZER: funções que salvam em arquivos de texto os parâmetros da malha e da 
// solução

int main(){
  /* MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  char filename_msh_x[] = "../P02/results_unmod/eom_x_iter_0000009781.dat";
  char filename_msh_y[] = "../P02/results_unmod/eom_y_iter_0000009781.dat";

  /* FLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  // sim_prmtrs config_fp;
  // config_fp.Ntype = 2;
  // config_fp.w = 2.;
  // config_fp.r = 1.;
  // config_fp.alpha_seq = 1;
  // config_fp.alpha = 1e-4;
  // config_fp.alpha_H = 1e-1;
  // config_fp.set_alpha_H = 0;
  // config_fp.M = 5;
  // config_fp.max_iter = 500000;
  // config_fp.qtimes = 72605;
  // config_fp.save_i_c = 0;
  // config_fp.save_last_only = 1;
  // config_fp.eps = 1.e-6;
  // char output_file_msh[] = "results/mesh";

  /* FLOW (solution) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  

  // print_2d_array_to_file(msh.JMAX,msh.IMAX,x,"./results/test_x.dat",0);
  // print_2d_array_to_file(msh.JMAX,msh.IMAX,y,"./results/test_y.dat",0);

  // free(msh.af_prmtrs);
  // free(msh.end_prmtrs);
  // free(config_msh.casename);
  // free(c_prmtrs.al);
  // free(c_prmtrs.bm);
  // free(c_prmtrs.cl);
  // free(c_prmtrs.dm);
  // free(c_prmtrs.ksi_l);
  // free(c_prmtrs.ksi_m);
  // free(c_prmtrs.eta_l);
  // free(c_prmtrs.eta_m);
  // free(x);
  // free(y);

  return 0;
}