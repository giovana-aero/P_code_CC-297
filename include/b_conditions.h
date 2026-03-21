#ifndef _lib_b_conditions_
#define _lib_b_conditions_

/*
Defined in the context of laplace2d
-gigiaero, 20/03/2026, 1436 hours
*/
typedef struct B_conditions_2d{
  char type;
  double val;
  int axis; // 1 for x, 2 for y
  int position; // index for the axis opposite to the one above
  int range[2];
}b_conditions_2d;

void apply_b_cs(int m, int n,double A[m][n],int num_b_cs,b_conditions_2d *b_c);
void dirichlet_rectangular_constant(int m,int n,double A[m][n],
                                    b_conditions_2d *b_c);

#endif