#include<stdio.h>
#include<stdlib.h>
#include"../include/bi_air_lib.h"
#include"../include/ecm_lib.h"
#include"../include/eom_lib.h"
#include"../include/mesh.h"

int main(){
  // Solution configurations
  sim_prmtrs config;
  config.Ntype = 2;
  config.w = 2.;
  config.r = 1.;
  config.alpha_seq = 1;
  config.alpha = 1e-4;
  config.alpha_H = 1e-1;
  config.set_alpha_H = 0;
  config.M = 5;
  config.max_iter = 500000;
  config.qtimes = 72605;
  config.save_i_c = 1;
  config.save_last_only = 1;
  config.eps = 1.e-6;
  char output_file[] = "results/eom";

  // Mesh parameters
  msh_prmtrs msh;
  msh.end_prmtrs = malloc(sizeof(double)*4);
  /* IMAX */
  // msh.IMAX = 93;
  msh.IMAX = 241;
  /* JMAX */
  // msh.JMAX = 15;
  msh.JMAX = 81;
  /* c */
  msh.c = 1.;
  /* end_prmtrs */
  msh.end_prmtrs[0] = 6.5*msh.c;
  msh.end_prmtrs[1] = 6.5*msh.c;
  msh.end_prmtrs[2] = msh.c*.5;
  // msh.end_prmtrs[2] = msh.c*.25;
  msh.end_prmtrs[3] = 0.;
  /* init_type */
  msh.init_type = 3;
  int init_only = 0; // Initialize only, do not solve elliptical mesh (1 or 0)
  int mesh_type = 1; // 1: eom; 2: ecm (ADI only)

  /* af_type */
  int n = 10; // cst - bernstein polynomial order
  msh.af_prmtrs = malloc(sizeof(double)*((n+2)*2 + 1));
  msh.af_type = 1;
  /* af_prmtrs (bi_air) */
  msh.af_prmtrs[0] = 0.1;
  /* af_prmtrs (naca4) */
  // msh.af_prmtrs[0] = 0.;
  // msh.af_prmtrs[1] = 0.;
  // msh.af_prmtrs[2] = 10.;
  // /* af_prmtrs (cst) */
  // msh.af_prmtrs[0] = (double) n;
  // msh.cst_foil = 3; 
  // cst_prmtrs(&msh);

  // P & Q control functions
  control_prmtrs c_prmtrs;
  int control_type = 0; // see set_control_prmtrs for list of configurations
  set_control_prmtrs(control_type,&c_prmtrs,&msh);

  // Solve
  config.casename = malloc(sizeof(char)*200);
  sprintf(config.casename,"%s",output_file);
  switch(mesh_type){
    case 1:
      evaluate_delta_form_eom(&config,&msh,&c_prmtrs,init_only);
      break;
    
    case 2:
      evaluate_delta_form_ecm(&config,&msh,&c_prmtrs,init_only);
      break;
    
    default:
      puts("invalid mesh_type");
      return 1;
  }
  
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

  return 0;
}