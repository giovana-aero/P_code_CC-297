#include"./eom_lib.h"

#ifndef _lib_pom_
#define _lib_pom_

void calc_xy_eta(double *xy_eta,int m,int n,double x[m][n],double y[m][n],
                 int i,int j,double S_eta);
void parabolic_mesh(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);

#endif