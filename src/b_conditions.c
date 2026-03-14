#include<stdio.h>
#include"../include/b_conditions.h"

typedef enum {vertical=0, horizontal=1} orientation;

/*
- gigiaero, 13/03/2026, 2223 hours
*/
void dirichlet_rectangular_constant(int m,int n,double A[m][n],double val,
                                    int *i_range,int fixed_i,int orientation){

  switch(orientation){
    case 0: // Horizontal
      for(int i=i_range[0];i<i_range[1];i++)
        A[fixed_i][i] = val;
      break;

    case 1: // Vertical
      for(int i=i_range[0];i<i_range[1];i++)
        A[i][fixed_i] = val;
      break;

    default: puts("Wrong orientation value");
  }
}