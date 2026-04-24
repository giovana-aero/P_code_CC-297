#include<stdio.h>
#include<stdlib.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/bi_air_lib.h"
#include"../include/eom_lib.h"
#include"../include/mesh.h"

int main(){

  double c = 1.; // Chord length

  msh_prmtrs msh;
  msh.af_prmtrs = malloc(sizeof(double));
  msh.end_prmtrs = malloc(sizeof(double)*4);
  // n
  msh.n = 93;
  // m
  msh.m = 15;
  // af_type
  msh.af_type = 1;
  // af_prmtrs
  msh.af_prmtrs[0] = 0.1;
  // end_prmtrs
  msh.end_prmtrs[0] = 6.5*c;
  msh.end_prmtrs[1] = 6.5*c;
  msh.end_prmtrs[2] = c*.5;
  msh.end_prmtrs[3] = 0.;
  // init_type
  msh.init_type = 1;
  


  double (*x)[msh.n] = calloc(msh.m,sizeof *x);
  double (*y)[msh.n] = calloc(msh.m,sizeof *y);

  double *x_tmp = malloc(sizeof(double)*msh.n);
  double *y_tmp = malloc(sizeof(double)*msh.n);

  ellipse(x_tmp,y_tmp,msh.end_prmtrs,msh.n);
  bi_air_shape(msh.af_prmtrs[0],)

  // save_mesh(msh.n,msh.n,x,y,"test");

  free(msh.af_prmtrs);
  free(msh.end_prmtrs);

  free(x);
  free(y);
  free(x_tmp);
  free(y_tmp);

  return 0;
}