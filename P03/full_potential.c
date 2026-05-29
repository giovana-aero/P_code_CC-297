#include<stdio.h>
#include<stdlib.h>
#include"../include/2d_arrays.h"
#include"../include/eom_lib.h"
#include"../include/mesh.h"

// nome padrão dos arquivos de malha: fullp_mesh_x.dat e fullp_mesh_y.dat
// init_only = 1 -> renomear a malha inicial pro nome padrão
// init_only = 0 -> renomear a malha final pro nome padrão

// A FAZER: funções que salvam em arquivos de texto os parâmetros da malha e da 
// solução

int main(){
  /* MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  int mesh_only = 1; // resolve mesh only, do not solve flow (1 or 0)
  int load_mesh = 1; // load mesh from file (1 or 0)
  char filename_msh_x[] = "../P02/results_unmod/eom_x_iter_0000009781.dat";
  char filename_msh_y[] = "../P02/results_unmod/eom_y_iter_0000009781.dat";

  // Solution configurations
  sim_prmtrs config_msh;
  config_msh.Ntype = 2;
  config_msh.w = 2.;
  config_msh.r = 1.;
  config_msh.alpha_seq = 1;
  config_msh.alpha = 1e-4;
  config_msh.alpha_H = 1e-1;
  config_msh.set_alpha_H = 0;
  config_msh.M = 5;
  config_msh.max_iter = 500000;
  config_msh.eps = 1.e-6;
  char output_file_msh[] = "results/fullp_mesh";

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
  int init_only = 0; // Initialize only, do not solve elliptical mesh (1 or 0)

  /* af_type */
  int n = 10; // cst - bernstein polynomial order
  msh.af_prmtrs = malloc(sizeof(double)*((n+2)*2 + 1));
  msh.af_type = 2;
  /* af_prmtrs (bi_air) */
  // msh.af_prmtrs[0] = 0.1;
  /* af_prmtrs (naca4) */
  msh.af_prmtrs[0] = 2.;
  msh.af_prmtrs[1] = 4.;
  msh.af_prmtrs[2] = 12.;
  /* af_prmtrs (cst) */
  // msh.af_prmtrs[0] = n;
  // int cst_foil = 1; // see cst_prmtrs for list of airfoils
  // cst_prmtrs(cst_foil,msh.af_prmtrs);

  // P & Q control functions
  control_prmtrs c_prmtrs;
  int control_type = 0; // see set_control_prmtrs for list of configurations

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

  /* MESH (solution) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  double (*x)[msh.IMAX] = calloc(msh.JMAX,sizeof *x);
  double (*y)[msh.IMAX] = calloc(msh.JMAX,sizeof *y);
  config_msh.casename = malloc(sizeof(char)*200);
  if(load_mesh){
    read_2d_array_from_file(msh.JMAX,msh.IMAX,x,filename_msh_x);
    read_2d_array_from_file(msh.JMAX,msh.IMAX,y,filename_msh_y);
  }
  else{
    set_control_prmtrs(control_type,&c_prmtrs,&msh);
    sprintf(config_msh.casename,"%s",output_file_msh);
    evaluate_delta_form_eom_fp(msh.JMAX,msh.IMAX,x,y,&config_msh,&msh,&c_prmtrs,
                               init_only);
  }

  /* FLOW (solution) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  if(!mesh_only){

  }

  print_2d_array_to_file(msh.JMAX,msh.IMAX,x,"./results/test_x.dat",0);
  print_2d_array_to_file(msh.JMAX,msh.IMAX,y,"./results/test_y.dat",0);

  free(msh.af_prmtrs);
  free(msh.end_prmtrs);
  free(config_msh.casename);
  free(c_prmtrs.al);
  free(c_prmtrs.bm);
  free(c_prmtrs.cl);
  free(c_prmtrs.dm);
  free(c_prmtrs.ksi_l);
  free(c_prmtrs.ksi_m);
  free(c_prmtrs.eta_l);
  free(c_prmtrs.eta_m);
  free(x);
  free(y);

  return 0;
}