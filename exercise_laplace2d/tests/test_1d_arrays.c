#include<stdio.h>
#include"../../include/1d_arrays.h"

int main(){

  int n = 7;

  double v[n];

  ones_1d(n,v);
  print_1d_array(n,v);
  putchar('\n');
  
  zeros_1d(n,v);
  print_1d_array(n,v);
  putchar('\n');

  fill_1d_array(n,v);
  print_1d_array(n,v);
  putchar('\n');

  return 0;
  
}