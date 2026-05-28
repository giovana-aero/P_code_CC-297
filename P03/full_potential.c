#include<stdio.h>
#include<stdlib.h>
#include"../include/eom_lib.h"
#include"../include/mesh.h"

int main(){
  int mesh_only = 1; // solve mesh only, do not solve flow (1 or 0)
  int load_mesh = 0; // load from file (1 or 0)

  /* MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
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
  config_msh.qtimes = 72605;
  config_msh.save_i_c = 1;
  config_msh.save_last_only = 1;
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
  // int mesh_type = 1; // 1: eom; 2: ecm (ADI only)

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
  if(load_mesh){

  }
  else{
    set_control_prmtrs(control_type,&c_prmtrs,&msh);
    config_msh.casename = malloc(sizeof(char)*200);
    sprintf(config_msh.casename,"%s",output_file_msh);
    // switch(mesh_type){
    //   case 1:
    evaluate_delta_form_eom(&config_msh,&msh,&c_prmtrs,init_only);
    //     break;
      
    //   case 2:
    //     evaluate_delta_form_ecm(&config,&msh,&c_prmtrs,init_only);
    //     break;
      
    //   default:
    //     puts("invalid mesh_type");
    //     return 1;
    // }
  }

  /* FLOW (solution) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  if(!mesh_only){

  }

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

  return 0;
}