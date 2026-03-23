#include<stdio.h>
#include<stdlib.h>
#include"../include/1d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/num_methods.h"

/*
- gigiaero, 20/03/2026, 2300 hours
*/
void apply_b_c(int m, int n,double A[m][n],int num_b_c,b_conditions_2d *b_c,
               double *x,double *y){

  for(int i=0;i<num_b_c;i++){
    if(b_c[i].axis != 1 && b_c[i].axis != 2){
      printf("apply_b_cs: Invalid axis specified at b_c[%d]\n",i);
      exit(28);
    }

    switch(b_c[i].type){
      case 'D': // Dirichlet, constant
        dirichlet_rectangular_constant(m,n,A,&b_c[i]);
        break;
      
      case 'd': // Dirichlet, variable
        dirichlet_rectangular_variable(m,n,A,&b_c[i]);
        break;
      
      case 'N': // Neumann (implemented with central differences)
        neumann_rectangular_constant(m,n,A,&b_c[i],x,y);
        break;

      default:
        printf("apply_b_c: invalid type at b_c[%d]\n",i);
        exit(28);
    }
  }
}

/*
- gigiaero, 13/03/2026, 2223 hours

modified to use b_conditions_2d structs
- gigiaero, 20/03/2026, 1439

loop range altered: now it's defined by <= instead of <
- gigiaero, 23/03/2026, 0945 hours

struct b_conditions_2d altered: now, for type 'D' boundary conditions, val must
be defined as a length 1 array
- gigiaero, 23/03/2026, 1007 hours
*/
void dirichlet_rectangular_constant(int m,int n,double A[m][n],
                                    b_conditions_2d *b_c){
  switch(b_c->axis){
    case 1: // Horizontal
      for(int i=b_c->range[0];i<=b_c->range[1];i++)
        A[b_c->position][i] = b_c->val[0];
      break;

    case 2: // Vertical
      for(int i=b_c->range[0];i<=b_c->range[1];i++)
        A[i][b_c->position] = b_c->val[0];
      break;
  }
}

/*
- gigiaero, 23/03/2026, 1048 hours
*/
void dirichlet_rectangular_variable(int m,int n,double A[m][n],
                                    b_conditions_2d *b_c){
  switch(b_c->axis){
    case 1: // Horizontal
      for(int i=b_c->range[0],idx=0;i<=b_c->range[1];i++,idx++)
        A[b_c->position][i] = b_c->val[idx];
      break;

    case 2: // Vertical
      for(int i=b_c->range[0],idx=0;i<=b_c->range[1];i++,idx++)
        A[i][b_c->position] = b_c->val[idx];
      break;
  }
}

void neumann_rectangular_constant(int m,int n,double A[m][n],
                                  b_conditions_2d *b_c,double *x, double *y){

  double (*tmp_A)[3] = calloc(3,sizeof *tmp_A);
  double *tmp_xy = malloc(sizeof(double)*3);
  double N,L_A;
  int old_A_rc_size;

  if(b_c->axis == 1)
    old_A_rc_size = n;
  else
    old_A_rc_size = m;

  double *old_A_rc = malloc(sizeof(double)*old_A_rc_size);

  switch(b_c->axis){
    case 1: // Horizontal
      if(b_c->position == 0){ // Down
        copy_1d_array_range(0,n,A[0],old_A_rc);

        tmp_xy[0] = -y[1];
        tmp_xy[1] = y[0];
        tmp_xy[2] = y[1];

        for(int i=b_c->range[0];i<=b_c->range[1];i++){
          tmp_A[1][0] = A[0][i-1];
          tmp_A[1][1] = A[0][i];
          tmp_A[1][2] = A[0][i+1];
          tmp_A[2][1] = A[1][i];
          tmp_A[0][1] = A[1][i] + 2.*delta_xy(tmp_xy,1)*b_c->val[0];
          // Reminder: the signal of the flux changes due to the normal vector
          
          scheme_der2_o2_central_var_deltas_xy(&L_A,3,3,tmp_A,x,tmp_xy,1,1);
          N_p_jacobi(&N,x,tmp_xy,i,1);
          A[0][i] = -L_A/N + old_A_rc[i];
        }
      }
      else if(b_c->position == m-1){  // Up
        // copy_1d_array_range(0,n,A[0],old_A_rc);
        for(int i=0;i<n;i++)
          old_A_rc[i] = A[b_c->position][i];
        
        // print_1d_array(n,old_A_rc);
        // for(int i=0;i<n;i++)
        //   printf("%f ",A[b_c->position][i]);
        // putchar('\n'); putchar('\n');

        tmp_xy[0] = y[m-2];
        tmp_xy[1] = y[m-1];
        tmp_xy[2] = y[m-1] + (y[m-1] - y[m-2]);

        // print_1d_array(3,tmp_xy);

        for(int i=b_c->range[0];i<=b_c->range[1];i++){
          tmp_A[1][0] = A[m-1][i-1];
          tmp_A[1][1] = A[m-1][i];
          tmp_A[1][2] = A[m-1][i+1];
          tmp_A[0][1] = A[m-2][i];
          tmp_A[2][1] = A[m-2][i] + 2.*delta_xy(tmp_xy,1)*b_c->val[0];

          // print_2d_array(m,n,A);
          // print_2d_array(3,3,tmp_A);
          // putchar('\n');
          
          scheme_der2_o2_central_var_deltas_xy(&L_A,3,3,tmp_A,x,tmp_xy,1,1);
          N_p_jacobi(&N,x,tmp_xy,i,1);
          A[m-1][i] = -L_A/N + old_A_rc[i];
        }
      }
      else{
        puts("Boundary condition placed in internal domain");
        exit(28);
      }
      
      break;

    case 2: // Vertical
      puts("b_conditions: case not yet implemented!");
      exit(28);
      // for(int j=0;j<m;j++)
      //   old_A_rc[j] = A[]
      // copy_1d_array_range(0,n,A[0],old_A_rc);

      // tmp_xy[0] = -x[1];
      // tmp_xy[1] = x[0];
      // tmp_xy[2] = x[1];

      // for(int i=b_c->range[0];i<=b_c->range[1];i++){
      //   tmp_A[1][0] = A[0][i-1];
      //   tmp_A[1][1] = A[0][i];
      //   tmp_A[1][2] = A[0][i+1];
      //   tmp_A[2][1] = A[1][i];
      //   tmp_A[0][1] = A[1][i] - 2.*delta_xy(tmp_xy,1)*b_c->val[0];
        
      //   scheme_der2_o2_central_var_deltas_xy(&L_A,3,3,tmp_A,x,tmp_xy,1,1);
      //   N_p_jacobi(&N,x,tmp_xy,i,1);
      //   A[0][i] = -L_A/N + old_A_rc[i];
      // }
      break;
  }

  free(tmp_A);
  free(tmp_xy);
  free(old_A_rc);

}