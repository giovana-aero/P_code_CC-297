#include"./eom_lib.h"

#ifndef _lib_pom_
#define _lib_pom_

void calc_A_line(int m,int n,double *A,double x[m][n],double y[m][n],int j);
void calc_B_line(int m,int n,double *A,double x[m][n],double y[m][n],int j);
void calc_C_line(int m,int n,double *A,double x[m][n],double y[m][n],int j);
void calc_xy_eta(double *xy_eta,int m,int n,double x[m][n],double y[m][n],
                 int i,int j,double S_eta);
void exspace(double *x,double a,int n);
void loc_ref_grid(int m,int n,double x[m][n],double y[m][n],double *eps_switch,
                  double *s,int j);
void parabolic_mesh(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);

#endif