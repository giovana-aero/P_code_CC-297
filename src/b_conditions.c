#include<stdio.h>
#include<stdlib.h>
#include"../include/b_conditions.h"

/*
-gigiaero, 20/03/2026, 2300 hours
*/
void apply_b_cs(int m, int n,double A[m][n],int num_b_cs,b_conditions_2d *b_c){
  for(int i=0;i<num_b_cs;i++){
    if(b_c[i].axis != 1 && b_c[i].axis != 2){
      printf("apply_b_cs: Invalid axis specified at b_c[%d]\n",i);
      exit(28);
    }

    switch(b_c[i].type){
      case 'D': // Dirichlet
        dirichlet_rectangular_constant(m,n,A,&b_c[i]);
        break;
      
      case 'N': // Neumann
        // pending!!!!1
        break;

      default:
        printf("apply_b_cs: invalid type at b_c[%d]\n",i);
        exit(28);
    }
  }
}

/*
-gigiaero, 13/03/2026, 2223 hours

modified to use b_conditions_2d structs
-gigiaero, 20/03/2026, 1439
*/
void dirichlet_rectangular_constant(int m,int n,double A[m][n],b_conditions_2d *b_c){

  switch(b_c->axis){
    case 1: // Horizontal
      for(int i=b_c->range[0];i<b_c->range[1];i++)
        A[b_c->position][i] = b_c->val;
      break;

    case 2: // Vertical
      for(int i=b_c->range[0];i<b_c->range[1];i++)
        A[i][b_c->position] = b_c->val;
      break;

    // default: puts("dirichlet_rectangular_constant: Wrong orientation value");
  }
}