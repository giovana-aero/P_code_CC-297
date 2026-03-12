#include<stdio.h>
#include"../../include/gigi_lib.h"

int main(){

  int n = 7;

  double v[n];

  // ones_1d(n,v);
  print_1d_array(n,v);
  putchar('\n');
  
  // zeros_1d(n,v);
  print_1d_array(n,v);
  putchar('\n');

  // fill_2d_array(m,n,A);
  // print_2d_array(m,n,A);
  // putchar('\n');

  return 0;
  
}