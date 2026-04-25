#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/bi_air_lib.h"
#include"../include/eom_lib.h"
#include"../include/mesh.h"

int main(){
  msh_prmtrs msh;
  msh.af_prmtrs = malloc(sizeof(double)*3);
  msh.end_prmtrs = malloc(sizeof(double)*4);
  // n
  msh.IMAX = 93;
  // m
  // msh.JMAX = 15;
  msh.JMAX = 30;
  // c
  msh.c = 1.;
  // af_type
  msh.af_type = 2;
  // af_prmtrs (bi_air)
  // msh.af_prmtrs[0] = 0.1;
  // af_prmtrs (naca4)
  msh.af_prmtrs[0] = 8.;
  msh.af_prmtrs[1] = 4.;
  msh.af_prmtrs[2] = 12;
  // end_prmtrs
  msh.end_prmtrs[0] = 6.5*msh.c;
  msh.end_prmtrs[1] = 6.5*msh.c;
  msh.end_prmtrs[2] = msh.c*.5;
  msh.end_prmtrs[3] = 0.;
  // init_type
  msh.init_type = 1;



  double (*x)[msh.IMAX] = calloc(msh.JMAX,sizeof *x);
  double (*y)[msh.IMAX] = calloc(msh.JMAX,sizeof *y);

  initialize_mesh(msh.JMAX,msh.IMAX,x,y,&msh);



  // double *v = malloc(sizeof(double)*4);
  // copy_1d_array_range(0,3,msh.end_prmtrs,v);

  // double ri = msh.end_prmtrs[0]/(msh.JMAX);
  // for(int j=msh.JMAX-1;j>=0;j--){
  //   ellipse(x[j],y[j],v,msh.IMAX,1);
  //   v[0] -= ri;
  //   v[1] -= ri;
  // }

  print_2d_array_to_file(msh.JMAX,msh.IMAX,x,"mesh_x.dat",0);
  print_2d_array_to_file(msh.JMAX,msh.IMAX,y,"mesh_y.dat",0);

  // double *x_tmp = malloc(sizeof(double)*msh.IMAX);
  // double *y_tmp = malloc(sizeof(double)*msh.IMAX);
  // ellipse(x_tmp,y_tmp,msh.end_prmtrs,msh.IMAX,1);

  // ellipse(x[1],y[1],msh.end_prmtrs,msh.IMAX);
  // bi_air_shape(msh.IMAX,msh.af_prmtrs[0],)

  // save_mesh(msh.IMAX,msh.IMAX,x_tmp,y_tmp,"test");

  free(msh.af_prmtrs);
  free(msh.end_prmtrs);

  free(x);
  free(y);
  // free(x_tmp);
  // free(y_tmp);

  return 0;
}