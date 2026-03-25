#ifndef _lib_b_conditions_
#define _lib_b_conditions_

/*
Defined in the context of laplace2d
- gigiaero, 20/03/2026, 1436 hours
*/
typedef struct B_conditions_2d{
  char type;
  double *val;
  int axis; // 1 for x, 2 for y
  int position; // index for the axis opposite to the one above
  int range[2];
}b_conditions_2d;

void apply_b_c(int m, int n,double A[m][n],int num_b_cs,b_conditions_2d *b_c,
               double *x,double *y);
void build_tmp_A_neumann_y_down(int m,int n,double A[m][n],double tmp_A[3][3],
                                double *tmp_y,double val,int i);
void dirichlet_rectangular_constant(int m,int n,double A[m][n],
                                    b_conditions_2d *b_c);
void dirichlet_rectangular_variable(int m,int n,double A[m][n],
                                    b_conditions_2d *b_c);
void neumann_rectangular_constant(int m,int n,double A[m][n],
                                  b_conditions_2d *b_c,double *x, double *y);

#endif