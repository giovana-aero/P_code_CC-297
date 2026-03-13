#include<stdio.h>
// #include"../../include/gigi_lib.so"
#include"../../include/2d_arrays.h"

int main(){

  int m = 5; // lines
  int n = 7; // columns

  double A[m][n];

  ones_2d(m,n,A);
  print_2d_array(m,n,A);
  putchar('\n');
  
  zeros_2d(m,n,A);
  print_2d_array(m,n,A);
  putchar('\n');

  fill_2d_array(m,n,A);
  print_2d_array(m,n,A);
  putchar('\n');

  return 0;
  
}