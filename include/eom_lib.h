#ifndef _lib_eom_
#define _lib_eom_

typedef struct mesh_parameters{
  int IMAX; // number of points (circunference)
  int JMAX; // number of points (normal)
  double c; // chord length
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

void calc_A(int m,int n,double A[m][n],double x[m][n],double y[m][n]);
void calc_B(int m,int n,double B[m][n],double x[m][n],double y[m][n]);
void calc_C(int m,int n,double C[m][n],double x[m][n],double y[m][n]);
void calc_D(int m,int n,double D[m][n],double x[m][n],double y[m][n]);
void cosspace(double *x,double xi,double xf,int n,int half);
void cst_airfoil(int n_pts,double *x,double *yu,double *yl,double *prmtrs,
                 double c);
void cst_prmtrs(int af,double *prmtrs);
void ellipse(double *x,double *y,double *prmtrs,int n,int invert_th);
void evaluate_delta_form_eom(sim_prmtrs *config,msh_prmtrs *msh,int init_only);
void init_af_bi_air(double *x,double *y,double *x_axis,int chord_n,
                    msh_prmtrs *msh);
void init_af_cst(double *x,double *y,double *x_axis,int chord_n,
                 msh_prmtrs *msh);
void init_af_naca4(double *x,double *y,double *x_axis,int chord_n,
                   msh_prmtrs *msh);
void init_type1(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void init_type2(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void init_type3(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void initialize_mesh(int m,int n,double x[m][n],double y[m][n],msh_prmtrs *msh);
void L_phi_eom(int m,int n,double L_phi_xy[m][n],double xy[m][n],double A[m][n],
               double B[m][n],double C[m][n],double D[m][n]);
void linspace(double *x,double xi,double xf,int n);
double max_thickness(int m,int n,double y[m][n]);
void naca4(int n,double *x,double *xu,double *xl,double *yu,double *yl,
           double *prmtrs);
void solve_adi_2d_rectangular_eom(int m,int n,double x[m][n],double y[m][n],
                                  sim_prmtrs *config);

#endif