#include<stdio.h>
#include"../include/1d_arrays.h"
#include"../include/2d_arrays.h"
#include"../include/b_conditions.h"
#include"../include/mesh.h"

int main(){

  // Mesh
  int m = 5;
  int n = 9;

  // Physical properties
  double alpha = 1.; // CONFIRMAR ISTO
  double lx = 2.;
  double ly = 2.;

  // Boundary conditions
  double Tu = 1.;
  double Td = 2.;
  double Tl = 3.;
  double Tr = 4.;
  int range_h[] = {1,n-1};
  int range_v[] = {1,m-1};

  // Initialize mesh and temperature array
  double x[n], y[m], T[m][n];
  double delta_x = lx/(n - 1.);
  double delta_y = ly/(m - 1.);
  zeros_2d(m,n,T);
  uniform_rectangular_mesh(m,n,delta_x,delta_y,x,y);

  // Initialize boundary conditions
  dirichlet_rectangular_constant(m,n,T,Tu,range_h,m-1,0);
  dirichlet_rectangular_constant(m,n,T,Td,range_h,0,0);
  dirichlet_rectangular_constant(m,n,T,Tl,range_v,0,1);
  dirichlet_rectangular_constant(m,n,T,Tr,range_v,n-1,1);

  // Solve


  // Print results to file
  

  // print_2d_array(m,n,T);

  return 0;
}