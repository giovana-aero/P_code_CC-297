#include<stdio.h>
#include<stdlib.h>
#include"../include/fullp_lib.h"

int main(){
  /* MESH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /* Biconvex airfoil, 10% thickness */
  // // Original mesh from P02, eom
  // char filename_msh_x[] = "../P02/eom_m15_n93_biair_t01/eom_x_iter_0000000047.dat";
  // char filename_msh_y[] = "../P02/eom_m15_n93_biair_t01/eom_y_iter_0000000047.dat";
  // int m = 15;
  // int n = 93;
  // // Arbitrarily refined, eom
  // char filename_msh_x[] = "../P02/eom_m60_n187_biair_t01/eom_x_iter_0000000069.dat";
  // char filename_msh_y[] = "../P02/eom_m60_n187_biair_t01/eom_y_iter_0000000069.dat";
  // int m = 60;
  // int n = 187;
  // Arbitrarily refined, with outer circle centered at 0.25c, eom
  char filename_msh_x[] = "../P02/eom_m60_n187_biair_t01_25c/eom_x_iter_0000000074.dat";
  char filename_msh_y[] = "../P02/eom_m60_n187_biair_t01_25c/eom_y_iter_0000000074.dat";
  int m = 60;
  int n = 187;
  // // Arbitrarily refined, eom
  // char filename_msh_x[] = "../P02/eom_m81_n241_biair_t01/eom_x_iter_0000000065.dat";
  // char filename_msh_y[] = "../P02/eom_m81_n241_biair_t01/eom_y_iter_0000000065.dat";
  // int m = 81;
  // int n = 241;
  // // Mesh 1, eom
  // char filename_msh_x[] = "../P02/eom_m11_n31_biair_t01/eom_x_iter_0000000168.dat";
  // char filename_msh_y[] = "../P02/eom_m11_n31_biair_t01/eom_y_iter_0000000168.dat";
  // int m = 11;
  // int n = 31;
  // // Mesh 2, eom
  // char filename_msh_x[] = "../P02/eom_m21_n61_biair_t01/eom_x_iter_0000000047.dat";
  // char filename_msh_y[] = "../P02/eom_m21_n61_biair_t01/eom_y_iter_0000000047.dat";
  // int m = 21;
  // int n = 61;
  // // Mesh 3, eom
  // char filename_msh_x[] = "../P02/eom_m41_n121_biair_t01/eom_x_iter_0000000064.dat";
  // char filename_msh_y[] = "../P02/eom_m41_n121_biair_t01/eom_y_iter_0000000064.dat";
  // int m = 41;
  // int n = 121;
  // // Mesh 3, eom, double outer radius
  // char filename_msh_x[] = "../P02/eom_m41_n121_biair_t01_outer_double/eom_x_iter_0000000112.dat";
  // char filename_msh_y[] = "../P02/eom_m41_n121_biair_t01_outer_double/eom_y_iter_0000000112.dat";
  // int m = 41;
  // int n = 121;
  // // Mesh 4, eom
  // char filename_msh_x[] = "../P02/eom_m81_n241_biair_t01/eom_x_iter_0000000065.dat";
  // char filename_msh_y[] = "../P02/eom_m81_n241_biair_t01/eom_y_iter_0000000065.dat";
  // int m = 81;
  // int n = 241;
  // // Mesh 4, pom
  // char filename_msh_x[] = "../P02/eom_m81_n241_biair_t01/pom_x_iter_0000000065.dat";
  // char filename_msh_y[] = "../P02/eom_m81_n241_biair_t01/pom_y_iter_0000000065.dat";
  // int m = 81;
  // int n = 241;

  /* NACA 0010 */
  // // Arbitrarily refined, eom
  // char filename_msh_x[] = "../P02/eom_m60_n187_naca0010/eom_x_iter_0000000069.dat";
  // char filename_msh_y[] = "../P02/eom_m60_n187_naca0010/eom_y_iter_0000000069.dat";
  // int m = 60;
  // int n = 187;
  // Mesh 3, eom
  // char filename_msh_x[] = "../P02/eom_m41_n121_naca0010/eom_x_iter_0000000064.dat";
  // char filename_msh_y[] = "../P02/eom_m41_n121_naca0010/eom_y_iter_0000000064.dat";
  // int m = 41;
  // int n = 121;

  /* NACA 0012 */
  // // Mesh 3, eom
  // char filename_msh_x[] = "../P02/eom_m41_n121_naca0012/eom_x_iter_0000000063.dat";
  // char filename_msh_y[] = "../P02/eom_m41_n121_naca0012/eom_y_iter_0000000063.dat";
  // int m = 41;
  // int n = 121;
  // // Mesh 4, pom
  // char filename_msh_x[] = "../P02/pom_m81_n241_naca0012/eom_x_initial.dat";
  // char filename_msh_y[] = "../P02/pom_m81_n241_naca0012/eom_y_initial.dat";
  // int m = 81;
  // int n = 241;

  /* NACA 64A010*/
  // // Mesh 3, eom
  // char filename_msh_x[] = "../P02/eom_m41_n121_naca64a010/eom_x_iter_0000000064.dat";
  // char filename_msh_y[] = "../P02/eom_m41_n121_naca64a010/eom_y_iter_0000000064.dat";
  // int m = 41;
  // int n = 121;

  /* supercritical TM-X-1831 */
  // // Mesh 3, eom
  // char filename_msh_x[] = "../P02/eom_m41_n121_supercritical_tm-x-1831/eom_x_iter_0000000107.dat";
  // char filename_msh_y[] = "../P02/eom_m41_n121_supercritical_tm-x-1831/eom_y_iter_0000000107.dat";
  // int m = 41;
  // int n = 121;
  // // Mesh 4, eom
  // char filename_msh_x[] = "../P02/eom_m81_n241_supercritical_tm-x-1831/eom_x_iter_0000000209.dat";
  // char filename_msh_y[] = "../P02/eom_m81_n241_supercritical_tm-x-1831/eom_y_iter_0000000209.dat";
  // int m = 81;
  // int n = 241;
  // // Mesh 4, pom
  // char filename_msh_x[] = "../P02/pom_m81_n241_supercritical_tm-x-1831/eom_x_initial.dat";
  // char filename_msh_y[] = "../P02/pom_m81_n241_supercritical_tm-x-1831/eom_y_initial.dat";
  // int m = 81;
  // int n = 241;

  /* FLOW ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  sim_prmtrs config;
  config.Ntype = 3;
  config.w = 1.;
  config.r = 1.;
  config.alpha_seq = 1;
  config.alpha = 1e-1;
  config.alpha_H = 2e-1;
  // config.alpha = 1e-1;
  // config.alpha_H = 2e-1;
  config.set_alpha_H = 0;
  config.M = 5;
  config.max_iter = 1000;
  config.qtimes = 1;
  config.save_i_c = 1;
  config.save_last_only = 1;
  config.eps = 1.e-6;
  // config.eps = 1.e-10;
  char output_file[] = "results/fullp";

  fullp_prmtrs fp_prmtrs;
  fp_prmtrs.alpha = 0.;
  fp_prmtrs.Ma = .84;
  fp_prmtrs.C = 1.1;
  fp_prmtrs.beta_sub = .3;
  fp_prmtrs.beta_super = 4.5;
  fp_prmtrs.lift = 1; // (not working)
  fp_prmtrs.Rg = .5;

  /* FLOW (solution) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  config.casename = malloc(sizeof(char)*200);
  sprintf(config.casename,"%s",output_file);
  evaluate_delta_form_fullp(m,n,&config,&fp_prmtrs,filename_msh_x,
                            filename_msh_y);

  free(config.casename);

  return 0;
}