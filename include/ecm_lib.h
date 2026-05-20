#include"./num_methods.h"
#include"./eom_lib.h"

#ifndef _lib_ecm_
#define _lib_ecm_

void evaluate_delta_form_ecm(sim_prmtrs *config,msh_prmtrs *msh,
                             control_prmtrs *c_prmtrs,int init_only);
void initialize_mesh_ecm(int m,int n,double x[m][n],double y[m][n],
                       msh_prmtrs *msh);
void half_ellipse(double *x,double *y,double *prmtrs,int n,int invert_th);

#endif