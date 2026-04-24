#ifndef _lib_eom_
#define _lib_eom_

typedef struct mesh_parameters{
  int n; // number of points (circunference)
  int m; // number of points (normal)
  /* 1-biconvex, 2-naca4, 3-cst */
  int af_type; 
  /*
  1-[t]
  2-[m,p,t]
  3-[(depends on the order of the bernstein polynomial)]
  */
  double *af_prmtrs;
  double *end_prmtrs; // [rx,ry,dx,dy]
  /*
  1-uniform (ie, always ellipses)
  2-linspace
  3-cosspace_half
  4-parabollic
  */
  int init_type;
}msh_prmtrs;

void cosspace(double *x,double xi,double xf,int n);
void ellipse(double *x,double *y,double *prmtrs,int n);
void linspace(double *x,double xi,double xf,int n);

#endif